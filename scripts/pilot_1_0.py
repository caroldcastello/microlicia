#!/usr/bin/env python3
"""Pilot 1.0: MCC tree, stamen morphology, and occurrence summaries.

Run from this directory or elsewhere:
  PILOT_INPUT_DIR=/path/to/input python3 pilot_1_0.py

Required files in PILOT_INPUT_DIR:
  teste_datacao_MCC.tree, an occurrence CSV, estames.csv
  occurrence input priority: ocorrencias.xlsx.csv, then ocorrencias.csv
"""
from __future__ import annotations
import csv, hashlib, json, math, os, re
import sys, platform, xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT_DIR=Path(os.environ.get('PILOT_INPUT_DIR','data/input'))
OUT=Path(os.environ.get('PILOT_OUTPUT_DIR','results'))
OUT.mkdir(parents=True,exist_ok=True)
TREE_PATH=INPUT_DIR/'teste_datacao_MCC.tree'
occurrence_env=os.environ.get('PILOT_OCCURRENCE_FILE')
if occurrence_env:
    OCC_PATH=Path(occurrence_env)
else:
    occurrence_candidates=[INPUT_DIR/'ocorrencias.xlsx.csv',INPUT_DIR/'ocorrencias.csv']
    OCC_PATH=next((p for p in occurrence_candidates if p.exists()),occurrence_candidates[0])
TRAIT_PATH=INPUT_DIR/'estames.csv'
XML_PATH=INPUT_DIR/'teste_datacao(1).xml'
DATING_REF=Path(os.environ.get('PILOT_DATING_REFERENCE','data/metadata/fonte_datação.json'))
FILTER_AUDIT=Path(os.environ.get('PILOT_FILTER_AUDIT','data/metadata/taxa_ocorrencias_excluidos_sem_correspondencia.csv'))

# MCC tip labels define the analysis reference. Candidate occurrence spellings
# are mapped to exact tree labels for matching only. This is not a claim of
# taxonomic synonymy; all original labels remain in the crosswalk and inputs.
OCCURRENCE_TO_TREE_TIP_ALIASES={
    'm_amblyssepala':'m_amblysepala',
    'm_amplexiaculis':'m_amplexicaulis',
    'm_aracnoidea':'m_arachnoidea',
    'm_bacchaoides':'m_baccharoides',
    'm_belinelloi':'m_bellineloi',
    'm_diamantinesis':'m_diamantinensis',
    'm_helvoa':'m_helvola',
    'm_hilairei':'m_hilarei',
}
GENUS_PREFIXES={
    'm':'Microlicia',
    'mar':'Marcetia',
    'p':'Poteranthera',
    'r':'Rhynchanthera',
    'rhe':'Rhexia',
    'si':'Siphanthera',
}
FIGURE_EPITHET_CORRECTIONS={
    # The submitted label has a typographic error; retain the original in the
    # input MCC and display the botanical spelling in the figure.
    'r_brachyrryncha':'Rhynchanthera brachyrhyncha',
}

# Tree parsing supports BEAST/NEXUS comments and a single Newick tree statement.
@dataclass
class Node:
    name:str|None=None
    children:list['Node']=field(default_factory=list)
    length:float=0.0
    ident:int=-1
    depth:float=0.0
    age:float=0.0

def parse_tree(path:Path)->Node:
    txt=path.read_text(encoding='utf-8',errors='replace')
    tr=re.search(r'(?is)\bTranslate\s+(.*?);',txt)
    translate={}
    if tr:
        for key,value in re.findall(r"(?m)^\s*(\d+)\s+([^,\s]+)\s*,?\s*$",tr.group(1)):
            translate[key]=value.strip("'\"")
    m=re.search(r'(?im)^\s*tree\s+[^=]+=(.*?);',txt,re.S)
    if not m: raise ValueError('Could not find a single NEXUS tree statement')
    s=m.group(1); pos=0
    def skip():
        nonlocal pos
        while pos<len(s):
            if s[pos].isspace(): pos+=1
            elif s[pos]=='[':
                end=s.find(']',pos+1); pos=len(s) if end<0 else end+1
            else: break
    def readtoken(stops):
        nonlocal pos
        start=pos
        while pos<len(s) and s[pos] not in stops: pos+=1
        return s[start:pos].strip()
    def parse():
        nonlocal pos
        skip(); children=[]; name=None
        if pos>=len(s): raise ValueError('Unexpected end of Newick tree')
        if s[pos]=='(':
            pos+=1
            while True:
                children.append(parse()); skip()
                if pos<len(s) and s[pos]==',': pos+=1; continue
                if pos<len(s) and s[pos]==')': pos+=1; break
                raise ValueError(f'Unexpected token near tree position {pos}')
            skip(); name=readtoken(':,();[') or None
        else:
            name=readtoken(':,();[') or None
        skip(); length=0.0
        if pos<len(s) and s[pos]==':':
            pos+=1; skip(); number=readtoken(',);[')
            length=float(number)
        return Node(name,children,length)
    root=parse()
    def resolve_translate(n):
        if not n.children and n.name in translate:n.name=translate[n.name]
        for child in n.children:resolve_translate(child)
    resolve_translate(root)
    nodes=[]; tips=[]
    def index(n,depth):
        n.ident=len(nodes); n.depth=depth; nodes.append(n)
        if not n.children: tips.append(n)
        for c in n.children: index(c,depth+c.length)
    index(root,0.0)
    root_age=max(t.depth for t in tips)
    for n in nodes:n.age=root_age-n.depth
    if any(not math.isfinite(n.length) or n.length<0 for n in nodes): raise ValueError('Invalid branch length')
    if len({t.name for t in tips})!=len(tips): raise ValueError('Duplicate tree tip labels')
    return root

def descendants(n):
    if not n.children:return [n.name]
    return sum((descendants(c) for c in n.children),[])

def read_occurrences(path):
    # Exports may be UTF-8 or Windows-1252; try UTF-8 first and fall back only
    # when the byte stream is not valid UTF-8.
    try:
        d=pd.read_csv(path,sep=';',encoding='utf-8-sig',dtype={'Taxon_original':'string'},keep_default_na=False)
    except UnicodeDecodeError:
        d=pd.read_csv(path,sep=';',encoding='cp1252',dtype={'Taxon_original':'string'},keep_default_na=False)
    needed={'Taxon_original','Latitude_decimal','Longitude_decimal'}
    if not needed.issubset(d.columns):raise ValueError(f'Occurrence columns missing: {needed-set(d.columns)}')
    d['taxon_raw']=d['Taxon_original'].astype(str).str.strip()
    d['taxon_normalized']=d['taxon_raw'].map(normalize_taxon_label)
    if 'taxon_tree' in d.columns:
        d['taxon_normalized']=d['taxon_tree'].astype(str).map(normalize_taxon_label)
    d['latitude']=pd.to_numeric(d['Latitude_decimal'],errors='coerce')
    d['longitude']=pd.to_numeric(d['Longitude_decimal'],errors='coerce')
    d['coord_valid']=d.latitude.between(-90,90)&d.longitude.between(-180,180)
    return d

def normalize_taxon_label(value):
    """Normalize formatting and a full Microlicia genus name to MCC-style labels."""
    x=str(value).strip().lower()
    x=re.sub(r'^microlicia[ _]+','m_',x)
    x=re.sub(r'^m[ _]+','m_',x)
    x=re.sub(r'\s+','_',x)
    x=re.sub(r'_+','_',x)
    return x

def read_traits(path):
    d=pd.read_csv(path,sep=';',encoding='utf-8-sig',dtype=str,keep_default_na=False)
    # The first field has an empty header and contains the taxon labels.
    first=d.columns[0]
    d=d.rename(columns={first:'taxon'})
    d['taxon']=d.taxon.astype(str).str.strip()
    d=d[d.taxon.ne('')].copy()
    expected=[
        ('trait1','dimorfismo estaminal','(0) dimorphic; (1) subisomorphic; (2) isomorphic'),
        ('trait2','padrão de coloração dos estames','(0) concolours; (1) bicolours'),
        ('trait3','esporângios das anteras','(0) tetrasporangiate; (1) polysporangiate')]
    if len(d.columns)<4:raise ValueError('Expected taxon plus three categorical traits')
    mapping={}
    for i,(key,name,description) in enumerate(expected,start=1):
        col=d.columns[i]; d=d.rename(columns={col:key}); mapping[key]={'name':name,'description':description}
        d[key]=pd.to_numeric(d[key].replace('',np.nan),errors='coerce')
    return d,mapping

def discrete_er(root, tip_states, states):
    k=len(states); idx={s:i for i,s in enumerate(states)}
    observed={t:tip_states.get(t) for t in tip_states}
    nodes=[]
    def collect(n):
        nodes.append(n)
        for c in n.children:collect(c)
    collect(root)
    def likelihood(logq,save=False):
        q=math.exp(float(logq)); Q=np.full((k,k),q); np.fill_diagonal(Q,-(k-1)*q)
        pmat={n.ident:expm(Q*n.length) for n in nodes if n is not root}
        L={}; logscale={}; msgs={}
        for n in reversed(nodes):
            if not n.children:
                state=observed.get(n.name)
                vec=np.ones(k) if state is None or (isinstance(state,float) and np.isnan(state)) else np.eye(k)[idx[state]]
                sc=0.0
            else:
                vec=np.ones(k); sc=0.0
                for c in n.children:
                    msg=pmat[c.ident]@L[c.ident]; msgs[c.ident]=msg
                    vec*=msg; sc+=logscale[c.ident]
            z=float(np.max(vec))
            if z<=0 or not math.isfinite(z):return (-math.inf,None) if save else -math.inf
            vec=vec/z; sc+=math.log(z); L[n.ident]=vec; logscale[n.ident]=sc
        root_prior=np.full(k,1/k)
        logl=math.log(float(root_prior@L[root.ident]))+logscale[root.ident]
        if not save:return logl
        posterior={root.ident:(root_prior*L[root.ident])/(root_prior@L[root.ident])}
        outside={root.ident:root_prior}
        for n in nodes:
            if not n.children:continue
            for c in n.children:
                context=outside[n.ident].copy()
                for sib in n.children:
                    if sib.ident!=c.ident: context*=msgs[sib.ident]
                out=context@pmat[c.ident]
                out=out/np.sum(out) if np.sum(out)>0 else np.full(k,1/k)
                outside[c.ident]=out
                w=out*L[c.ident]; posterior[c.ident]=w/w.sum()
        return logl, (q,Q,posterior)
    opt=minimize_scalar(lambda x:-likelihood(x),bounds=(-9,7),method='bounded',options={'xatol':1e-8})
    logl, details=likelihood(opt.x,save=True)
    q,Q,posterior=details
    node_rows=[]
    for n in nodes:
        if n.children:
            row={'node_id':n.ident+1,'is_root':n is root,'age_tree_units':n.age}
            for i,s in enumerate(states):row[f'P_state_{s}']=float(posterior[n.ident][i])
            row['most_supported_state']=states[int(np.argmax(posterior[n.ident]))]
            node_rows.append(row)
    root_prob=posterior[root.ident]
    return {'q_per_tree_unit':q,'log_likelihood':logl,'root_probs':root_prob,'node_rows':node_rows,'states':states}

def write_csv(df,path):df.to_csv(path,index=False,encoding='utf-8-sig')

def sha256(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def clone_subtree_tips(node, keep):
    """Return the induced subtree for `keep`, collapsing unary nodes."""
    if not node.children:
        if node.name not in keep:return None
        return Node(node.name,[],node.length,-1,node.depth,node.age)
    children=[clone_subtree_tips(c,keep) for c in node.children]
    children=[c for c in children if c is not None]
    if not children:return None
    if len(children)==1:
        child=children[0]
        child.length+=node.length
        return child
    return Node(node.name,children,node.length,-1,node.depth,node.age)

def reindex_subtree(node):
    i=0
    def walk(n):
        nonlocal i
        n.ident=i;i+=1
        for child in n.children:walk(child)
    walk(node)
    node.length=0.0
    return node

def display_species_name(taxon):
    if taxon in FIGURE_EPITHET_CORRECTIONS:return FIGURE_EPITHET_CORRECTIONS[taxon]
    if '_' not in taxon:return taxon.replace('_',' ')
    genus_code,epithet=taxon.split('_',1)
    genus=GENUS_PREFIXES.get(genus_code,genus_code)
    return f'{genus} {epithet}'

def plot_dated_tree(root, path_png, path_pdf):
    """Plot the dated MCC with every tip, using the original MCC topology."""
    nodes=[]
    def collect(n):
        nodes.append(n)
        for c in n.children:collect(c)
    collect(root)
    tips=[n for n in nodes if not n.children]
    labels={n.name:display_species_name(n.name) for n in tips}
    # Keep the tip order encoded in the MCC instead of sorting species names.
    # This preserves the tree's visual clade order and makes the branching easier
    # to follow from the root to the tips.
    ordered=tips
    y_by_id={n.ident:i for i,n in enumerate(ordered)}
    clade_tips={}
    def collect_clade_tips(n):
        if not n.children:result={n.name}
        else:
            result=set().union(*(collect_clade_tips(c) for c in n.children))
        clade_tips[n.ident]=result
        return result
    collect_clade_tips(root)
    def set_y(n):
        if not n.children:return y_by_id[n.ident]
        ys=[set_y(c) for c in n.children]
        y_by_id[n.ident]=sum(ys)/len(ys)
        return y_by_id[n.ident]
    set_y(root)
    root_depth=min(n.depth for n in nodes)
    root_depth=min(root_depth,root.depth)
    tree_height=max(n.depth for n in tips)-root.depth
    fig_h=max(16.0,0.23*len(tips))
    fig,ax=plt.subplots(figsize=(18,fig_h))
    for n in nodes:
        if not n.children:continue
        child_y=[y_by_id[c.ident] for c in n.children]
        ax.vlines(n.depth,min(child_y),max(child_y),color='#a0a8a3',lw=.55,zorder=1)
        for c in n.children:
            is_microlicia_clade=all(t.startswith('m_') for t in clade_tips[c.ident])
            edge_color='#2d7250' if is_microlicia_clade else '#65727c'
            ax.hlines(y_by_id[c.ident],n.depth,c.depth,color=edge_color,lw=.72,zorder=1)
    for tip in ordered:
        color='#2d7250' if tip.name.startswith('m_') else '#46515b'
        ax.text(tip.depth+tree_height*.008,y_by_id[tip.ident],labels[tip.name],va='center',ha='left',fontsize=7.2,fontstyle='italic',color=color)
    root_age=root.age
    depth_max=max(n.depth for n in tips)
    ticks=np.linspace(root.depth,depth_max,6)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{root_age-(x-root.depth):.1f}' for x in ticks])
    ax.set_xlabel('Idade antes do presente (unidades numéricas da MCC)')
    ax.set_yticks([])
    ax.set_title(f'Árvore datada MCC — {len(tips)} espécies; {len([t for t in tips if t.name.startswith("m_")])} de Microlicia',loc='left',fontsize=14,pad=14)
    ax.scatter([root.depth],[y_by_id[root.ident]],s=22,marker='o',color='#222222',zorder=3)
    ax.annotate('raiz',xy=(root.depth,y_by_id[root.ident]),xytext=(tree_height*.015,y_by_id[root.ident]-4),fontsize=8,color='#222222',arrowprops={'arrowstyle':'-','color':'#222222','lw':.6})
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0],[0],color='#2d7250',lw=2,label='Microlicia'),Line2D([0],[0],color='#46515b',lw=2,label='Outros táxons da árvore')],loc='lower right',frameon=False)
    ax.set_xlim(root.depth-tree_height*.015,depth_max+tree_height*.42)
    ax.set_ylim(-1,len(tips))
    ax.invert_yaxis()
    ax.grid(axis='x',alpha=.2,lw=.5)
    fig.tight_layout()
    fig.savefig(path_png,dpi=300,bbox_inches='tight')
    fig.savefig(path_pdf,bbox_inches='tight')
    plt.close(fig)

if '--tree-figure-only' in sys.argv:
    if not TREE_PATH.exists():raise FileNotFoundError(f'Missing dated MCC tree: {TREE_PATH}')
    tree_for_figure=parse_tree(TREE_PATH)
    plot_dated_tree(tree_for_figure,OUT/'dated_tree_all_taxa.png',OUT/'dated_tree_all_taxa.pdf')
    print(f'Wrote {OUT / "dated_tree_all_taxa.png"} and {OUT / "dated_tree_all_taxa.pdf"}')
    raise SystemExit(0)

# Load inputs and perform structural checks.
for p in (TREE_PATH,TRAIT_PATH):
    if not p.exists():raise FileNotFoundError(f'Missing input: {p}')
root_full=parse_tree(TREE_PATH)
if XML_PATH.exists():
    xml_root=ET.parse(XML_PATH).getroot()
    dating_units=sorted({n.attrib['units'] for n in xml_root.iter() if n.tag in {'yuleModel','constantSize'} and 'units' in n.attrib})
    calibrations=[]
    for prior in xml_root.iter('normalPrior'):
        statistic=prior.find('statistic')
        if statistic is not None:
            calibrations.append({'statistic':statistic.attrib.get('idref',''),'mean':float(prior.attrib['mean']),'stdev':float(prior.attrib['stdev'])})
else:
    dating_ref=json.loads(DATING_REF.read_text(encoding='utf-8'))
    dating_units=dating_ref['xml_declared_units_observed']
    calibrations=dating_ref['calibration_priors']
allnodes=[]
def nodes_walk(n):allnodes.append(n);[nodes_walk(c) for c in n.children]
nodes_walk(root_full)
tips=[n for n in allnodes if not n.children]; tipset={n.name for n in tips}
root_age=max(n.depth for n in tips); root_range=max(n.depth for n in tips)-min(n.depth for n in tips)
microlicia_tipset={t for t in tipset if t.startswith('m_')}
root=clone_subtree_tips(root_full,microlicia_tipset)
if root is None:raise ValueError('No Microlicia tips were found in the MCC tree')
root=reindex_subtree(root)
micro_nodes=[]
def micro_nodes_walk(n):micro_nodes.append(n);[micro_nodes_walk(c) for c in n.children]
micro_nodes_walk(root)
micro_tips=[n for n in micro_nodes if not n.children]
micro_tipset={n.name for n in micro_tips}
plot_dated_tree(root_full,OUT/'dated_tree_all_taxa.png',OUT/'dated_tree_all_taxa.pdf')
occ_input=read_occurrences(OCC_PATH)
occ_input['taxon_match_exact']=occ_input.taxon_normalized
occ_input['taxon_match']=occ_input.taxon_normalized.replace(OCCURRENCE_TO_TREE_TIP_ALIASES)
occ_input['name_rule']=np.where(
    occ_input.taxon_normalized.isin(OCCURRENCE_TO_TREE_TIP_ALIASES),
    'candidate spelling mapped to MCC tip; taxonomic validation pending',
    np.where(occ_input.taxon_normalized.isin(tipset),'exact match after formatting normalization','no MCC match; excluded from tree-based analyses')
)
occ_input['taxon_matched_full_tree']=occ_input.taxon_match.isin(tipset)
occ_input['taxon_matched_microlicia']=occ_input.taxon_match.isin(microlicia_tipset)
traits_all,trait_meta=read_traits(TRAIT_PATH)
traits_all['taxon_raw']=traits_all.taxon
traits_all['taxon']=traits_all.taxon.map(normalize_taxon_label)
traits_all['matched_mcc_tip']=traits_all.taxon.isin(tipset)
traits_all['included_microlicia']=traits_all.taxon.isin(microlicia_tipset)
if traits_all.taxon.duplicated().any():raise ValueError('Duplicate trait taxon rows; resolve before analysis')
traits=traits_all[traits_all.taxon.isin(microlicia_tipset)].copy()
if not traits.taxon.isin(microlicia_tipset).all():raise ValueError('Trait filtering to Microlicia failed')
occ=occ_input[occ_input.taxon_matched_microlicia].copy()
if not len(occ):raise ValueError('No occurrence records matched Microlicia tips in the MCC after spacing-only normalization')
trait_taxa=set(traits.taxon)

# Species crosswalk and unmatched occurrence taxa.
occ_taxa_all=occ_input.groupby(['taxon_raw','taxon_normalized','taxon_match','name_rule'],dropna=False).agg(n_records=('taxon_raw','size'),n_valid_coordinates=('coord_valid','sum')).reset_index()
occ_taxa_all['matched_full_tree']=occ_taxa_all.taxon_match.isin(tipset)
occ_taxa_all['included_in_microlicia_analyses']=occ_taxa_all.taxon_match.isin(microlicia_tipset)
write_csv(occ_taxa_all,OUT/'occurrence_taxon_crosswalk.csv')
write_csv(occ_taxa_all,OUT/'name_reconciliation_occurrences.csv')
trait_name_audit=traits_all[['taxon_raw','taxon','matched_mcc_tip','included_microlicia']].copy()
trait_name_audit['name_rule']=np.where(trait_name_audit.matched_mcc_tip,'exact match after formatting normalization','no MCC match; excluded from tree-based analyses')
write_csv(trait_name_audit,OUT/'name_reconciliation_stamen_traits.csv')
unmatched=occ_taxa_all[~occ_taxa_all.matched_full_tree].copy().sort_values(['n_records','taxon_raw'],ascending=[False,True])
write_csv(unmatched,OUT/'occurrence_taxa_unmatched.csv')
outgroup_occurrence_taxa=occ_taxa_all[occ_taxa_all.matched_full_tree & ~occ_taxa_all.included_in_microlicia_analyses].copy()
write_csv(outgroup_occurrence_taxa,OUT/'occurrence_outgroups_excluded_from_analyses.csv')
unmatched_audit=unmatched

# Coordinate and per-species summary, limited to names with an exact tip after spacing-only normalization.
valid=occ[occ.coord_valid].copy()
matched_valid=valid[valid.taxon_matched_microlicia].copy()
species_geo=matched_valid.groupby('taxon_match').agg(
    n_records=('taxon_match','size'),n_unique_coordinates=('latitude',lambda x: len(set(zip(x,matched_valid.loc[x.index,'longitude'])))),
    latitude_mean=('latitude','mean'),latitude_median=('latitude','median'),latitude_min=('latitude','min'),latitude_max=('latitude','max'),
    longitude_mean=('longitude','mean'),longitude_median=('longitude','median'),longitude_min=('longitude','min'),longitude_max=('longitude','max')).reset_index()
species_geo['latitude_span_degrees']=species_geo.latitude_max-species_geo.latitude_min
species_geo['longitude_span_degrees']=species_geo.longitude_max-species_geo.longitude_min
write_csv(species_geo,OUT/'occurrence_species_summary_tree_matched.csv')
write_csv(matched_valid,OUT/'occurrence_records_tree_matched_valid.csv')
write_csv(occ[~occ.coord_valid],OUT/'occurrence_records_without_valid_coordinate.csv')
tree_taxa_without_occ=sorted(microlicia_tipset-set(occ_taxa_all.loc[occ_taxa_all.included_in_microlicia_analyses,'taxon_match']))
write_csv(pd.DataFrame({'microlicia_taxon_without_matched_occurrence':tree_taxa_without_occ}),OUT/'microlicia_taxa_without_occurrences.csv')

# Trait frequencies and labels.
trait_labels={
 'trait1':{0:'dimórficos',1:'subisomórficos',2:'isomórficos'},
 'trait2':{0:'concolores',1:'bicolores'},
 'trait3':{0:'tetrasporangiadas',1:'poliesporangiadas'}}
trait_freq=[]; recon_results={}
for col in trait_meta:
    counts=traits[col].value_counts(dropna=True).sort_index()
    for state,count in counts.items():
        trait_freq.append({'trait':trait_meta[col]['name'],'trait_column':col,'state_code':int(state),'state_label':trait_labels[col].get(int(state),str(state)),'n_species':int(count),'proportion_of_nonmissing':float(count/counts.sum())})
    state_map=dict(zip(traits.taxon,traits[col]))
    states=sorted(int(x) for x in traits[col].dropna().unique())
    recon=discrete_er(root,state_map,states); recon_results[col]=recon
    rows=[]
    for row in recon['node_rows']:
        row['trait']=trait_meta[col]['name']; row['trait_column']=col
        row['state_labels']='; '.join(f"{s}={trait_labels[col].get(s,str(s))}" for s in states)
        rows.append(row)
    write_csv(pd.DataFrame(rows),OUT/f'ancestral_reconstruction_{col}.csv')
    for sidx,state in enumerate(states):
        trait_freq.append({'trait':trait_meta[col]['name'],'trait_column':col,'state_code':state,'state_label':trait_labels[col].get(state,str(state)),'n_species':int(counts.get(state,0)),'proportion_of_nonmissing':float(counts.get(state,0)/counts.sum()),'root_posterior':float(recon['root_probs'][sidx])})
write_csv(pd.DataFrame(trait_freq).drop_duplicates(['trait_column','state_code'],keep='last'),OUT/'stamen_trait_frequencies_and_root_posterior.csv')

# Joined species-level descriptive comparison of stamen categories and coordinate centroids.
geo_lookup=species_geo.set_index('taxon_match')
joined=traits.copy()
joined['in_tree']=joined.taxon.isin(micro_tipset)
for c in ['latitude_mean','longitude_mean','n_records','latitude_span_degrees','longitude_span_degrees']:
    joined[c]=joined.taxon.map(geo_lookup[c])
write_csv(joined,OUT/'stamen_traits_with_occurrence_summaries.csv')
comparison=[]
for col,meta in trait_meta.items():
    for state,label in trait_labels[col].items():
        sub=joined[(joined[col]==state)&joined.latitude_mean.notna()]
        for variable in ['latitude_mean','longitude_mean','n_records','latitude_span_degrees','longitude_span_degrees']:
            vals=sub[variable].dropna().astype(float)
            comparison.append({'trait':meta['name'],'trait_column':col,'state_code':state,'state_label':label,'variable':variable,'n_species':len(vals),'median':float(vals.median()) if len(vals) else np.nan,'mean':float(vals.mean()) if len(vals) else np.nan,'minimum':float(vals.min()) if len(vals) else np.nan,'maximum':float(vals.max()) if len(vals) else np.nan})
write_csv(pd.DataFrame(comparison),OUT/'exploratory_trait_geography_descriptives.csv')

# Figures.
plt.figure(figsize=(8,7))
plt.scatter(matched_valid.longitude,matched_valid.latitude,s=7,alpha=.28,color='#2f7d52',label='Microlicia')
plt.xlabel('Longitude (graus decimais; sem projeção aplicada)');plt.ylabel('Latitude (graus decimais)');plt.title('Registros de ocorrência de Microlicia com coordenadas válidas')
plt.legend(frameon=False);plt.grid(alpha=.18);plt.tight_layout();plt.savefig(OUT/'occurrence_coordinates_scatter.png',dpi=180);plt.close()
fig,axs=plt.subplots(1,3,figsize=(12,4))
for ax,(col,meta) in zip(axs,trait_meta.items()):
    c=traits[col].value_counts().sort_index(); labels=[trait_labels[col].get(int(s),str(s)) for s in c.index]
    ax.bar(labels,c.values,color=['#507f63','#93b49c','#d6b35d'][:len(c)])
    ax.set_title(meta['name'].replace('_',' '));ax.set_ylabel('Número de espécies');ax.tick_params(axis='x',rotation=20)
fig.tight_layout();fig.savefig(OUT/'stamen_trait_state_counts.png',dpi=180);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(11,4))
for ax,(col,meta) in zip(axs,list(trait_meta.items())[:2]):
    data=[joined.loc[(joined[col]==state)&joined.latitude_mean.notna(),'latitude_mean'].values for state in sorted(trait_labels[col])]
    names=[trait_labels[col][state] for state in sorted(trait_labels[col])]
    ax.boxplot(data,tick_labels=names,showfliers=False);ax.set_title(meta['name'].replace('_',' '));ax.set_ylabel('Latitude média por espécie (°)');ax.tick_params(axis='x',rotation=20)
fig.tight_layout();fig.savefig(OUT/'exploratory_trait_latitude.png',dpi=180);plt.close(fig)

# Run manifest and report.
coord_complete=int(occ.coord_valid.sum()); coord_bad=int(len(occ)-coord_complete)
occ_taxa_micro=occ_taxa_all[occ_taxa_all.included_in_microlicia_analyses]
matched_unique=int(occ_taxa_micro.taxon_match.nunique())
unmatched_unique=int((~occ_taxa_all.matched_full_tree).sum())
outgroup_unique=int(occ_taxa_all.loc[occ_taxa_all.matched_full_tree & ~occ_taxa_all.included_in_microlicia_analyses,'taxon_match'].nunique())
matched_record_count=len(occ)
matched_valid_record_count=len(matched_valid)
tree_check=pd.DataFrame([{
    'n_tips':len(tips),'n_unique_tips':len(tipset),'n_nodes':len(allnodes),
    'n_microlicia_tips':len(micro_tipset),'n_microlicia_analysis_tips':len(micro_tips),
    'root_to_tip_range_tree_units':root_range,'is_ultrametric_tolerance_1e-8':root_range<1e-8,
    'root_height_tree_units':root_age,'microlicia_crown_age_tree_units':root.age,'minimum_branch_length':min(n.length for n in allnodes if n is not root_full),
    'maximum_branch_length':max(n.length for n in allnodes if n is not root_full),
    'xml_model_units':';'.join(dating_units),
    'calibration_priors':'; '.join(f"{c['statistic']}={c['mean']}±{c['stdev']}" for c in calibrations)
}])
write_csv(tree_check,OUT/'tree_check.csv')
overview=pd.DataFrame([
    {'measure':'occurrence_records_input','value':len(occ_input)},
    {'measure':'occurrence_records_in_microlicia_analyses','value':len(occ)},
    {'measure':'occurrence_records_for_non_microlicia_tree_tips_excluded','value':int((occ_input.taxon_matched_full_tree & ~occ_input.taxon_matched_microlicia).sum())},
    {'measure':'occurrence_records_unmatched_to_full_tree','value':int((~occ_input.taxon_matched_full_tree).sum())},
    {'measure':'unique_occurrence_taxa_input','value':int(occ_input.taxon_raw.nunique())},
    {'measure':'microlicia_tips_without_occurrences','value':len(tree_taxa_without_occ)},
    {'measure':'coordinate_pairs_in_numeric_bounds','value':coord_complete},
    {'measure':'coordinate_pairs_missing_or_out_of_bounds','value':coord_bad},
    {'measure':'microlicia_occurrence_taxa_matched_to_tree_after_spacing_rule','value':matched_unique},
    {'measure':'occurrence_taxa_for_non_microlicia_tips_excluded','value':outgroup_unique},
    {'measure':'occurrence_taxa_unmatched_to_full_tree','value':unmatched_unique},
    {'measure':'records_matched_to_tree','value':matched_record_count},
    {'measure':'matched_records_with_coordinates_in_numeric_bounds','value':matched_valid_record_count},
    {'measure':'duplicate_coordinate_rows_vs_unique_coordinate_pairs','value':int(coord_complete-valid[['latitude','longitude']].drop_duplicates().shape[0])}
])
write_csv(overview,OUT/'occurrence_overview.csv')
recon_summary=[]
for col,meta in trait_meta.items():
    recon=recon_results[col]; states=recon['states']
    recon_summary.append({'trait':meta['name'],'states':states,'q_per_tree_unit':recon['q_per_tree_unit'],'log_likelihood':recon['log_likelihood'],'root_probabilities':{trait_labels[col].get(s,str(s)):float(recon['root_probs'][i]) for i,s in enumerate(states)}})
metadata={
    'pilot_version':'1.2',
    'run_time_utc':datetime.now(timezone.utc).isoformat(),
    'software':{'python':sys.version.split()[0],'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':__import__('scipy').__version__,'matplotlib':matplotlib.__version__},
    'inputs':{p.name:{'sha256':sha256(p),'bytes':p.stat().st_size} for p in (TREE_PATH,OCC_PATH,TRAIT_PATH,*( [XML_PATH] if XML_PATH.exists() else [DATING_REF] ),*( [FILTER_AUDIT] if FILTER_AUDIT.exists() else [] ))},
    'analysis_scope':'Microlicia only; all non-Microlicia tips are excluded from comparative analyses after dating',
    'tree':{'n_tips_all_taxa':len(tips),'n_unique_tips_all_taxa':len(tipset),'n_microlicia_tips':len(microlicia_tipset),'n_non_microlicia_tips':len(tipset)-len(microlicia_tipset),'root_age_tree_units':root_age,'microlicia_crown_age_tree_units':root.age,'root_to_tip_range':root_range,'xml_declared_model_units':dating_units,'calibration_priors':calibrations,'absolute_scale_interpretation':'Ma intended by matching calibration values to Reginato et al. 2022 Table 3 (Ma); XML label says years and must be reconciled before final inference','tree_tip_labels_used_as_analysis_reference':True,'tree_tip_labels_modified':False},
    'occurrence':{'n_records_input':len(occ_input),'n_records_microlicia_analyses':len(occ),'n_records_non_microlicia_tree_tips_excluded':int((occ_input.taxon_matched_full_tree & ~occ_input.taxon_matched_microlicia).sum()),'n_records_unmatched_to_full_tree':int((~occ_input.taxon_matched_full_tree).sum()),'n_coordinate_pairs_valid_in_microlicia':coord_complete,'n_coordinate_pairs_needing_review_in_microlicia':coord_bad,'n_unique_taxa_in_input':int(occ_input.taxon_raw.nunique()),'n_microlicia_taxa_matched':matched_unique,'n_unmatched_taxa_to_full_tree':unmatched_unique,'n_microlicia_tips_without_occurrences':len(tree_taxa_without_occ),'n_matched_records':matched_record_count,'n_matched_records_with_valid_coordinates':matched_valid_record_count,'duplicate_coordinate_pairs_global_within_microlicia':int(coord_complete-valid[['latitude','longitude']].drop_duplicates().shape[0]),'occurrence_to_tree_candidate_aliases':OCCURRENCE_TO_TREE_TIP_ALIASES,'alias_mapping_is_taxonomic_synonymy_claim':False,'name_audit_output':'name_reconciliation_occurrences.csv'},
    'stamen':{'n_rows_input_with_taxon':len(traits_all),'n_microlicia_rows':len(traits),'n_microlicia_taxa':int(traits.taxon.nunique()),'n_exact_matches_to_microlicia_tree':int(traits.taxon.isin(micro_tipset).sum()),'all_three_traits_nonmissing_for_microlicia_rows':bool(traits[list(trait_meta)].notna().all(axis=None))},
    'ancestral_reconstruction_model':'Equal-rates CTMC (ER), maximum likelihood; MCC tree only; uniform root prior',
    'ancestral_reconstruction':recon_summary,
    'figures':{'dated_tree_all_taxa_png':'dated_tree_all_taxa.png','dated_tree_all_taxa_pdf':'dated_tree_all_taxa.pdf','tip_count':len(tips),'labels':'Full genus and species names; Microlicia highlighted separately'},
    'limitations':['Only species belonging to Microlicia are included in occurrence and stamen analyses; tips from other genera remain in the full dated-tree figure.','Unmatched input taxa and occurrences for non-Microlicia tips are excluded from the Microlicia analyses and listed in separate crosswalk outputs.','Coordinate records marked invalid or missing in the supplied standardized CSV are excluded from spatial summaries; raw names and coordinates remain in the inputs.','Geographic output is descriptive in decimal degrees without map projection, area calculation, environmental rasters, or region coding.','This run uses the MCC tree only; uncertainty across posterior trees is not propagated.','The intended scale is provisionally interpreted as Ma because the calibration values match Reginato et al. (2022), Table 3; the XML unit label says years and should be reconciled with the executed BEAST settings.','Exploratory summaries by stamen state are descriptive and do not establish adaptation or causation.']
}
(OUT/'run_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')

lines=[
'# Microlicia — Piloto de análises 1.2','',
'## Dados e execução','',
'As análises de ocorrências, estames e reconstrução ancestral incluem somente espécies de *Microlicia*. Uma figura separada apresenta a MCC completa, com os 153 táxons rotulados. Os nomes da árvore MCC definem a referência. Os rótulos originais de ocorrência e estames são preservados em tabelas de auditoria; os oito mapeamentos de grafia candidatos às pontas da MCC são usados apenas para cruzamento, sem afirmar sinonímia taxonômica e sem sobrescrever entradas.','',
'## Resultados de validação','',
f'- Árvore completa: {len(tips)} pontas únicas, sendo {len(microlicia_tipset)} espécies de *Microlicia* e {len(tips)-len(microlicia_tipset)} pontas de outros gêneros; ultramétrica (amplitude raiz–ponta = {root_range:.3g}); altura da raiz = {root_age:.6f} unidades numéricas da árvore. O XML disponível declara units={dating_units}; médias/priores avaliados: {", ".join(str(c["mean"]) for c in calibrations)}.',
'- [Inferência] A escala pretendida é Ma: 29,47 (Marcetieae), 30,03 (Rhexieae) e 22,39 (Lavoisiereae) correspondem às idades da Tabela 3 do capítulo de Reginato et al. (2022), cuja coluna está em Ma. O XML, porém, declara `years`; essa divergência deve ser reconciliada com a configuração efetivamente executada no BEAST antes de publicação ou reexecução da datação.',
f'- Estames: {len(traits)} táxons de *Microlicia* após o filtro; {int(traits.taxon.isin(micro_tipset).sum())} correspondem às pontas da subárvore de *Microlicia*.',
f'- Ocorrências: {len(occ_input)} registros na entrada; {len(occ)} registros de *Microlicia* foram incluídos; {int((occ_input.taxon_matched_full_tree & ~occ_input.taxon_matched_microlicia).sum())} registros de outros gêneros com pontas na MCC foram excluídos das análises; {unmatched_unique} táxons não corresponderam à árvore completa.',
f'- Entre os registros incluídos de *Microlicia*, {coord_complete} pares de coordenadas estão completos e dentro dos limites numéricos de latitude/longitude; {coord_bad} precisam de revisão.',
f'- Cobertura ligada à subárvore de *Microlicia*: {matched_record_count} registros de {matched_unique} táxons; {matched_valid_record_count} registros têm coordenadas válidas.',
'', '## Estames e reconstrução ancestral','',
'As frequências observadas estão em `stamen_trait_frequencies_and_root_posterior.csv`. A reconstrução ancestral foi feita na subárvore de *Microlicia*, por máxima verossimilhança sob um modelo de transições iguais (ER), com probabilidades iguais para os estados na raiz. As probabilidades por nó estão nos três arquivos `ancestral_reconstruction_trait*.csv`.','',
*([f"- Probabilidade posterior estimada na raiz da subárvore de *Microlicia* para **{meta['name']}**: " + '; '.join(f"{state} = {prob:.3f}" for state,prob in recon['root_probabilities'].items()) + '.' for meta,recon in zip(trait_meta.values(),recon_summary)]),
'[Inferência] A interpretação das taxas por Ma depende de resolver a divergência entre a tabela de calibração e o rótulo `years` do XML. A reconstrução depende do modelo e não incorpora variação entre árvores.','',
'## Ocorrências e associação descritiva','',
'`occurrence_species_summary_tree_matched.csv` resume número de registros, coordenadas únicas, centroides aritméticos e amplitude em graus por espécie de *Microlicia*. `stamen_traits_with_occurrence_summaries.csv` cruza caracteres e resumos geográficos. `exploratory_trait_geography_descriptives.csv` apresenta estatísticas descritivas por estado morfológico, sem teste causal. O critério de coordenada válida é numérico (latitude entre −90 e 90, longitude entre −180 e 180).','',
'## Figuras','',
'- `dated_tree_all_taxa.png` e `dated_tree_all_taxa.pdf`: árvore datada MCC completa, com as 153 espécies identificadas por binômios; *Microlicia* aparece destacada em verde. O eixo usa as unidades numéricas da MCC enquanto a escala das calibrações é reconciliada.',
'- `occurrence_coordinates_scatter.png`: pontos de *Microlicia* em longitude/latitude decimal, sem projeção ou limites políticos.',
'- `stamen_trait_state_counts.png`: frequências dos estados dos três caracteres de *Microlicia*.',
'- `exploratory_trait_latitude.png`: distribuição descritiva da latitude média por estado nos dois primeiros caracteres de *Microlicia*.',
'', '## Conferência de nomes e limites do Piloto 1.2','',
'- `name_reconciliation_occurrences.csv` registra rótulo original, normalização de formato, rótulo da MCC usado no cruzamento, regra aplicada e inclusão/exclusão. `name_reconciliation_stamen_traits.csv` faz a conferência da matriz de estames. Nomes sem correspondência permanecem na auditoria e não entram nas análises ancoradas na árvore.',
'- As análises comparativas excluem pontas de gêneros diferentes de *Microlicia*; a figura da MCC mostra todas as pontas da árvore.',
'- Táxons da entrada sem correspondência com a árvore e registros de outros gêneros constam em tabelas de auditoria separadas; não entram nos resumos de *Microlicia*.',
'- Nove registros de *Microlicia* sem par de coordenadas numérico válido não entram nos resumos espaciais.',
'- Não há CRS documentado, áreas biogeográficas codificadas nem camadas ambientais neste piloto. A execução usa a MCC; a amostra posterior não foi processada nesta versão.',
'- A reconstrução ER é exploratória; modelos alternativos e sensibilidade precisam ser avaliados antes de inferência final.',
'', '## Arquivos de saída','',
'`run_metadata.json` registra hashes, versões e limites. O arquivo original da MCC e os arquivos de entrada não são sobrescritos.'
]
(OUT/'relatorio_piloto_1_2.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({'output_dir':str(OUT),'tips':len(tips),'root_age_tree_units':root_age,'root_to_tip_range':root_range,'occurrence_records':len(occ),'coordinates_valid':coord_complete,'coordinates_review':coord_bad,'occurrence_taxa':int(occ.taxon_raw.nunique()),'tree_matched_taxa':matched_unique,'unmatched_taxa':unmatched_unique,'matched_records':matched_record_count,'matched_valid_records':matched_valid_record_count,'trait_taxa':len(traits),'trait_exact_tree_matches':int(traits.taxon.isin(tipset).sum()),'trait_reconstructions':recon_summary},ensure_ascii=False,indent=2))
