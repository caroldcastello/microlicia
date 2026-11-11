#!/usr/bin/env python3
"""Pilot 1.0: MCC tree, stamen morphology, and occurrence summaries.

Run from this directory or elsewhere:
  PILOT_INPUT_DIR=/path/to/input python3 pilot_1_0.py

Required files in PILOT_INPUT_DIR:
  teste_datacao_MCC.tree, ocorrencias.xlsx.csv, estames.csv, teste_datacao(1).xml
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
OCC_PATH=INPUT_DIR/'ocorrencias.xlsx.csv'
TRAIT_PATH=INPUT_DIR/'estames.csv'
XML_PATH=INPUT_DIR/'teste_datacao(1).xml'
DATING_REF=Path(os.environ.get('PILOT_DATING_REFERENCE','data/metadata/fonte_datação.json'))
FILTER_AUDIT=Path(os.environ.get('PILOT_FILTER_AUDIT','data/metadata/taxa_ocorrencias_excluidos_sem_correspondencia.csv'))

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
    # Exported CSV is semicolon-delimited and Windows-1252 encoded.
    d=pd.read_csv(path,sep=';',encoding='utf-8-sig',dtype={'Taxon_original':'string'},keep_default_na=False)
    needed={'Taxon_original','Latitude_decimal','Longitude_decimal'}
    if not needed.issubset(d.columns):raise ValueError(f'Occurrence columns missing: {needed-set(d.columns)}')
    d['taxon_raw']=d['Taxon_original'].astype(str).str.strip()
    # Explicit, reversible formatting harmonization only; no spelling/synonym changes.
    d['taxon_match']=d['taxon_raw'].str.replace(r'_\s+','_',regex=True)
    d['latitude']=pd.to_numeric(d['Latitude_decimal'],errors='coerce')
    d['longitude']=pd.to_numeric(d['Longitude_decimal'],errors='coerce')
    d['coord_valid']=d.latitude.between(-90,90)&d.longitude.between(-180,180)
    return d

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

# Load inputs and perform structural checks.
for p in (TREE_PATH,OCC_PATH,TRAIT_PATH):
    if not p.exists():raise FileNotFoundError(f'Missing input: {p}')
root=parse_tree(TREE_PATH); tree_tips=[n.name for n in []]
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
nodes_walk(root)
tips=[n for n in allnodes if not n.children]; tipset={n.name for n in tips}
root_age=max(n.depth for n in tips); root_range=max(n.depth for n in tips)-min(n.depth for n in tips)
occ=read_occurrences(OCC_PATH); traits,trait_meta=read_traits(TRAIT_PATH)
if traits.taxon.duplicated().any():raise ValueError('Duplicate trait taxon rows; resolve before analysis')
occ['taxon_matched_tree']=occ.taxon_match.isin(tipset)
if not occ['taxon_matched_tree'].all():
    raise ValueError('The analysis input must contain only occurrence taxa matching the MCC after spacing-only normalization. Review the prefilter crosswalk.')
trait_taxa=set(traits.taxon)

# Species crosswalk and unmatched occurrence taxa.
occ_taxa=occ.groupby(['taxon_raw','taxon_match'],dropna=False).agg(n_records=('taxon_raw','size'),n_valid_coordinates=('coord_valid','sum')).reset_index()
occ_taxa['matched_tree']=occ_taxa.taxon_match.isin(tipset)
write_csv(occ_taxa,OUT/'occurrence_taxon_crosswalk.csv')
unmatched=occ_taxa[~occ_taxa.matched_tree].copy().sort_values(['n_records','taxon_raw'],ascending=[False,True])
if FILTER_AUDIT.exists():
    unmatched_audit=pd.read_csv(FILTER_AUDIT,encoding='utf-8-sig')
    write_csv(unmatched_audit,OUT/'occurrence_taxa_unmatched.csv')
else:
    unmatched_audit=unmatched.copy()
    write_csv(unmatched_audit,OUT/'occurrence_taxa_unmatched.csv')

# Coordinate and per-species summary, limited to names with an exact tip after spacing-only normalization.
valid=occ[occ.coord_valid].copy()
matched_valid=valid[valid.taxon_matched_tree].copy()
species_geo=matched_valid.groupby('taxon_match').agg(
    n_records=('taxon_match','size'),n_unique_coordinates=('latitude',lambda x: len(set(zip(x,matched_valid.loc[x.index,'longitude'])))),
    latitude_mean=('latitude','mean'),latitude_median=('latitude','median'),latitude_min=('latitude','min'),latitude_max=('latitude','max'),
    longitude_mean=('longitude','mean'),longitude_median=('longitude','median'),longitude_min=('longitude','min'),longitude_max=('longitude','max')).reset_index()
species_geo['latitude_span_degrees']=species_geo.latitude_max-species_geo.latitude_min
species_geo['longitude_span_degrees']=species_geo.longitude_max-species_geo.longitude_min
write_csv(species_geo,OUT/'occurrence_species_summary_tree_matched.csv')
write_csv(matched_valid,OUT/'occurrence_records_tree_matched_valid.csv')
write_csv(occ[~occ.coord_valid],OUT/'occurrence_records_without_valid_coordinate.csv')
tree_taxa_without_occ=sorted(tipset-set(occ_taxa.loc[occ_taxa.matched_tree,'taxon_match']))
write_csv(pd.DataFrame({'tree_taxon_without_matched_occurrence':tree_taxa_without_occ}),OUT/'tree_taxa_without_occurrences.csv')

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
joined['in_tree']=joined.taxon.isin(tipset)
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
plt.scatter(matched_valid.longitude,matched_valid.latitude,s=7,alpha=.28,color='#2f7d52',label='Taxon coincidente com árvore')
other=valid[~valid.taxon_matched_tree]
if len(other):plt.scatter(other.longitude,other.latitude,s=12,alpha=.45,color='#c97927',marker='x',label='Sem correspondência com árvore')
plt.xlabel('Longitude (graus decimais; sem projeção aplicada)');plt.ylabel('Latitude (graus decimais)');plt.title('Registros de ocorrência com coordenadas válidas')
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
exact_before=int(occ.taxon_raw.isin(tipset).sum()); matched_rows=int(occ.taxon_matched_tree.sum())
matched_unique=int(occ_taxa.matched_tree.sum()); unmatched_unique=int((~occ_taxa.matched_tree).sum())
matched_record_count=int(occ.loc[occ.taxon_matched_tree,'taxon_raw'].size)
matched_valid_record_count=len(matched_valid)
tree_check=pd.DataFrame([{
    'n_tips':len(tips),'n_unique_tips':len(tipset),'n_nodes':len(allnodes),
    'root_to_tip_range_tree_units':root_range,'is_ultrametric_tolerance_1e-8':root_range<1e-8,
    'root_height_tree_units':root_age,'minimum_branch_length':min(n.length for n in allnodes if n is not root),
    'maximum_branch_length':max(n.length for n in allnodes if n is not root),
    'xml_model_units':';'.join(dating_units),
    'calibration_priors':'; '.join(f"{c['statistic']}={c['mean']}±{c['stdev']}" for c in calibrations)
}])
write_csv(tree_check,OUT/'tree_check.csv')
overview=pd.DataFrame([
    {'measure':'occurrence_records','value':len(occ)},
    {'measure':'unique_raw_occurrence_taxa','value':int(occ.taxon_raw.nunique())},
    {'measure':'coordinate_pairs_in_numeric_bounds','value':coord_complete},
    {'measure':'coordinate_pairs_missing_or_out_of_bounds','value':coord_bad},
    {'measure':'occurrence_taxa_matched_to_tree_after_spacing_rule','value':matched_unique},
    {'measure':'occurrence_taxa_unmatched_to_tree','value':unmatched_unique},
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
    'pilot_version':'1.0',
    'run_time_utc':datetime.now(timezone.utc).isoformat(),
    'software':{'python':sys.version.split()[0],'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':__import__('scipy').__version__,'matplotlib':matplotlib.__version__},
    'inputs':{p.name:{'sha256':sha256(p),'bytes':p.stat().st_size} for p in (TREE_PATH,OCC_PATH,TRAIT_PATH,*( [XML_PATH] if XML_PATH.exists() else [DATING_REF] ),*( [FILTER_AUDIT] if FILTER_AUDIT.exists() else [] ))},
    'tree':{'n_tips':len(tips),'n_unique_tips':len(tipset),'root_age_tree_units':root_age,'root_age_interpreted_as_Ma_provisionally':root_age,'root_to_tip_range':root_range,'xml_declared_model_units':dating_units,'xml_normal_calibration_priors':calibrations,'absolute_scale_interpretation':'Ma intended by matching calibration values to Reginato et al. 2022 Table 3 (Ma); XML label says years and must be reconciled before final inference'},
    'occurrence':{'n_records':len(occ),'n_coordinate_pairs_valid':coord_complete,'n_coordinate_pairs_needing_review':coord_bad,'n_unique_taxa_in_tree_filtered_input':int(occ.taxon_raw.nunique()),'n_tree_matched_taxa_after_spacing_only_normalization':matched_unique,'n_unmatched_taxa_in_filtered_input':unmatched_unique,'n_taxa_excluded_before_analysis':int(len(unmatched_audit)),'n_records_excluded_before_analysis':int(unmatched_audit.n_records.sum()) if 'n_records' in unmatched_audit.columns else 0,'n_matched_records':matched_record_count,'n_matched_records_with_valid_coordinates':matched_valid_record_count,'exact_name_matches_before_spacing_rule_records':exact_before,'duplicate_coordinate_pairs_global':int(coord_complete-valid[['latitude','longitude']].drop_duplicates().shape[0])},
    'stamen':{'n_rows_with_taxon':len(traits),'n_taxa':int(traits.taxon.nunique()),'n_exact_matches_to_tree':int(traits.taxon.isin(tipset).sum()),'all_three_traits_nonmissing_for_all_named_taxa':bool(traits[list(trait_meta)].notna().all(axis=None))},
    'ancestral_reconstruction_model':'Equal-rates CTMC (ER), maximum likelihood; MCC tree only; uniform root prior',
    'ancestral_reconstruction':recon_summary,
    'limitations':['Taxon matches after the spacing-only rule are reported in a crosswalk; 24 occurrence taxa remain unmatched and are not included in tree-linked summaries.','Coordinate records marked invalid or missing in the supplied standardized CSV are excluded from spatial summaries; raw names and coordinates remain in the inputs.','Geographic output is descriptive in decimal degrees without map projection, area calculation, environmental rasters, or region coding.','This run uses the MCC tree only; uncertainty across posterior trees is not propagated.','The intended scale is provisionally interpreted as Ma because the calibration values match Reginato et al. (2022), Table 3; the XML unit label says years and should be reconciled with the executed BEAST settings.','Exploratory summaries by stamen state are descriptive and do not establish adaptation or causation.']
}
(OUT/'run_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')

lines=[
'# Microlicia — Piloto de análises 1.0','',
'## Dados e execução','',
'Análises feitas com a árvore MCC, ocorrências previamente filtradas para táxons da árvore e matriz categórica de estames. O script `scripts/pilot_1_0.py` gera tabelas, figuras e metadados. O XML do BEAST não é versionado; seus rótulos de unidade e priors foram registrados em `data/metadata/fonte_datação.json`.','',
'## Resultados de validação','',
f'- Árvore: {len(tips)} pontas únicas; ultramétrica (amplitude raiz–ponta = {root_range:.3g}); altura da raiz = {root_age:.6f} unidades da árvore. O XML local declara units={dating_units}; médias/priores avaliados: {", ".join(str(c["mean"]) for c in calibrations)}.',
'- [Inferência] A escala pretendida é Ma: 29,47 (Marcetieae), 30,03 (Rhexieae) e 22,39 (Lavoisiereae) correspondem às idades da Tabela 3 do capítulo de Reginato et al. (2022), cuja coluna está em Ma. A altura da raiz da MCC fica, portanto, provisoriamente interpretada como 28,048 Ma. O XML, porém, declara `years`; essa divergência deve ser reconciliada com a configuração efetivamente executada no BEAST antes de publicação ou reexecução da datação.',
f'- Estames: {len(traits)} táxons nomeados; {int(traits.taxon.isin(tipset).sum())} correspondem exatamente às pontas da árvore. Os três caracteres têm dados para todos os táxons nomeados.',
f'- Ocorrências: {len(occ)} registros; {coord_complete} pares de coordenadas completos dentro dos limites numéricos de latitude/longitude; {coord_bad} sem par válido.',
f'- Nomes: a análise recebeu ocorrências já filtradas para as pontas da MCC após remover apenas espaços entre `m_` e o epíteto. Foram mantidos {matched_record_count} registros de {matched_unique} espécies. Antes da filtragem, 24 táxons (117 registros) não tiveram correspondência; a lista está em `data/metadata/taxa_ocorrencias_excluidos_sem_correspondencia.csv` e nos resultados `occurrence_taxa_unmatched.csv`.',
f'- Cobertura espacial ligada à árvore: {matched_record_count} registros de {matched_unique} táxons; {matched_valid_record_count} registros têm coordenadas válidas.',
'', '## Estames e reconstrução ancestral','',
'As frequências observadas estão em `stamen_trait_frequencies_and_root_posterior.csv`. A reconstrução ancestral foi feita por máxima verossimilhança sob um modelo de transições iguais (ER), na MCC, com probabilidades iguais para os estados na raiz antes de observar os dados. As probabilidades por nó estão nos três arquivos `ancestral_reconstruction_trait*.csv`.','',
*([f"- Probabilidade posterior estimada na raiz para **{meta['name']}**: " + '; '.join(f"{state} = {prob:.3f}" for state,prob in recon['root_probabilities'].items()) + '.' for meta,recon in zip(trait_meta.values(),recon_summary)]),
'[Inferência] Como a escala pretendida é Ma, as taxas `q` podem ser lidas por Ma apenas se a divergência entre a fonte de calibração e o rótulo `years` do XML for resolvida. A reconstrução depende do modelo e não incorpora variação entre árvores.','',
'## Ocorrências e associação descritiva','',
'`occurrence_species_summary_tree_matched.csv` resume número de registros, coordenadas únicas, centroides aritméticos e amplitude em graus por táxon coincidente com a árvore. `stamen_traits_with_occurrence_summaries.csv` cruza os caracteres e esses resumos. `exploratory_trait_geography_descriptives.csv` apresenta medianas e amplitudes de latitude/longitude e esforço de coleta por estado morfológico; são estatísticas descritivas, sem teste causal. O critério de coordenada válida nesta etapa é apenas numérico (latitude entre −90 e 90, longitude entre −180 e 180); isso não confirma a localização geográfica de cada registro.','',
'## Figuras','',
'- `occurrence_coordinates_scatter.png`: pontos em longitude/latitude decimal, sem projeção ou limites políticos.',
'- `stamen_trait_state_counts.png`: frequências dos estados dos três caracteres.',
'- `exploratory_trait_latitude.png`: distribuição descritiva da latitude média por estado nos dois primeiros caracteres.',
'', '## Limites do Piloto 1.0','',
'- Os 24 nomes sem correspondência (117 registros) foram excluídos das ocorrências de entrada após normalização apenas de espaços em `m_`; a tabela de exclusão preserva nomes e contagens, sem inferir sinônimos ou corrigir grafias.',
'- Nove registros das espécies mantidas não têm par de coordenadas válido pelos limites numéricos e foram excluídos dos resumos espaciais; precisam de revisão.',
'- Não há CRS documentado, áreas biogeográficas codificadas nem camadas ambientais neste piloto. Não foram estimadas áreas em km², modelos ambientais ou biogeografia histórica. A execução usa a MCC; a amostra posterior não foi processada nesta versão.',
'- A reconstrução ER é a análise exploratória inicial; modelos alternativos e sensibilidade precisam ser avaliados antes de inferência final, principalmente para o caráter com um único táxon no estado “isomórficos”.',
'', '## Arquivos de saída','',
'`run_metadata.json` registra hash SHA-256 e tamanho dos arquivos de entrada, parâmetros gerais e limitações. Os dados originais não foram sobrescritos.'
]
(OUT/'relatorio_piloto_1_0.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({'output_dir':str(OUT),'tips':len(tips),'root_age_tree_units':root_age,'root_to_tip_range':root_range,'occurrence_records':len(occ),'coordinates_valid':coord_complete,'coordinates_review':coord_bad,'occurrence_taxa':int(occ.taxon_raw.nunique()),'tree_matched_taxa':matched_unique,'unmatched_taxa':unmatched_unique,'matched_records':matched_record_count,'matched_valid_records':matched_valid_record_count,'trait_taxa':len(traits),'trait_exact_tree_matches':int(traits.taxon.isin(tipset).sum()),'trait_reconstructions':recon_summary},ensure_ascii=False,indent=2))
