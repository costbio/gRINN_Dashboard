import os
import dash
from dash import Dash, dcc, html, dash_table, Input, Output, State, no_update
import pandas as pd
import plotly.graph_objects as go
import dash_molstar
from dash_molstar.utils import molstar_helper
from dash_molstar.utils.representations import Representation
import networkx as nx

# File paths
data_dir = os.path.join(os.path.dirname(__file__), 'test_data', 'prot_lig_1')
pdb_path = os.path.join(data_dir, 'system_dry.pdb')
total_csv = os.path.join(data_dir, 'energies_intEnTotal.csv')
traj_xtc = os.path.join(data_dir, 'traj_superposed.xtc')

# Load and transform interaction energy data
total_df = pd.read_csv(total_csv)
total_df['Pair'] = total_df['res1'] + '-' + total_df['res2']
cols2drop = [
    'Unnamed: 0','res1_index','res2_index','res1_chain','res2_chain',
    'res1_resnum','res2_resnum','res1_resname','res2_resname'
]
total_long = (
    total_df
    .drop(columns=cols2drop + ['res1', 'res2'])
    .melt(id_vars=['Pair'], var_name='Frame', value_name='Energy')
)
total_long['Energy'] = pd.to_numeric(total_long['Energy'], errors='coerce')
total_long = total_long[total_long['Energy'].notna()].copy()

# Determine frame range
df_frames = pd.to_numeric(total_long['Frame'], errors='coerce').dropna().astype(int)
frame_min, frame_max = int(df_frames.min()), int(df_frames.max())

# Residue list
first_res_list = sorted(total_df['res1'].unique())

# Molecular visualization setup
cartoon = Representation(type='cartoon', color='uniform')
cartoon.set_color_params({'value': 0xD3D3D3})
chainA = molstar_helper.get_targets(chain='A')
component = molstar_helper.create_component(label='Protein', targets=[chainA], representation=cartoon)
topo = molstar_helper.parse_molecule(pdb_path, component=component)
coords = molstar_helper.parse_coordinate(traj_xtc)

def get_full_trajectory():
    return molstar_helper.get_trajectory(topo, coords)

initial_traj = get_full_trajectory()

# Build graph helper
def build_graph(frame, include_cov, cutoff):
    df_f = total_long[total_long['Frame'].astype(int) == frame]
    G = nx.Graph()
    for res in first_res_list:
        G.add_node(res)
    for _, row in df_f.iterrows():
        r1, r2 = row['Pair'].split('-')
        e = row['Energy']
        if abs(e) >= cutoff:
            G.add_edge(r1, r2, weight=abs(e))
    if 'include' in include_cov:
        for i in range(len(first_res_list) - 1):
            G.add_edge(first_res_list[i], first_res_list[i+1], weight=0.0)
    return G

# App layout
app = Dash(__name__)
app.layout = html.Div([
    html.H1("gRINN Workflow Results", style={'textAlign': 'center'}),
    html.Div(style={'display': 'flex', 'height': '100vh', 'gap': '5px'}, children=[
        # Left Panel: Tabs
        html.Div(style={'width': '65%', 'border': '1px solid #CCC', 'padding': '10px', 'boxSizing': 'border-box'}, children=[
            dcc.Tabs(id='main-tabs', value='tab-pairwise', children=[
                # Pairwise Energies
                dcc.Tab(label='Pairwise Energies', value='tab-pairwise', children=[
                    html.Div(style={'display': 'flex', 'height': 'calc(100vh - 50px)', 'gap': '2px'}, children=[
                        html.Div(style={'minWidth': '160px', 'maxWidth': '200px', 'overflowY': 'auto'}, children=[
                            html.H4("Select First Residue"),
                            dash_table.DataTable(
                                id='first_residue_table',
                                columns=[{'name': 'Residue', 'id': 'Residue'}],
                                data=[{'Residue': r} for r in first_res_list],
                                row_selectable='single',
                                style_table={'height': 'calc(100vh - 200px)', 'overflowY': 'scroll'}
                            )
                        ]),
                        html.Div(style={'width': '200px', 'maxWidth': '200px', 'overflowY': 'auto'}, children=[
                            html.H4("Select Second Residue & IE"),
                            dash_table.DataTable(
                                id='second_residue_table',
                                columns=[{'name': 'Residue', 'id': 'Residue'},{'name': 'IE [kcal/mol]', 'id': 'IE'}],
                                data=[],
                                row_selectable='single',
                                style_table={'height': 'calc(100vh - 200px)', 'overflowY': 'scroll'}
                            )
                        ]),
                        html.Div(style={'flex': '2', 'paddingLeft': '10px'}, children=[
                            dcc.Graph(id='pair_energy_graph', style={'height': 'calc(100vh - 100px)'})
                        ])
                    ])
                ]),
                # Interaction Energy Matrix
                dcc.Tab(label='Interaction Energy Matrix', value='tab-matrix', children=[
                    dcc.Graph(id='matrix_heatmap', style={'height': 'calc(100vh - 100px)'})
                ]),
                # Network Analysis
                dcc.Tab(label='Network Analysis', value='tab-network', children=[
                    html.Div([
                        # Controls
                        html.Div(style={'display': 'flex', 'alignItems': 'center', 'paddingBottom': '10px'}, children=[
                            dcc.Checklist(
                                id='include_covalent_edges',
                                options=[{'label': 'Include covalent bonds as edges', 'value': 'include'}],
                                value=['include'],
                                style={'marginRight': '20px'}
                            ),
                            html.Label("Edge addition energy cutoff (kcal/mol): "),
                            dcc.Input(id='energy_cutoff', type='number', value=1.0, step=0.1, style={'width': '80px', 'marginRight': '20px'}),
                            html.Button('Update Network', id='update_network_btn', n_clicks=0),
                            html.Button('Export Network to File...', id='export_network_btn', n_clicks=0, style={'marginLeft': '10px'})
                        ]),
                        # Subtabs
                        dcc.Tabs(id='network-tabs', value='residue-metrics', children=[
                            dcc.Tab(label='Residue Metrics', value='residue-metrics', children=[
                                html.Div(style={'display': 'flex', 'gap': '20px', 'height': '60vh', 'overflowY': 'auto'}, children=[
                                    dcc.Graph(id='degree_centrality', style={'flex': '1'}),
                                    dcc.Graph(id='betweenness_centrality', style={'flex': '1'}),
                                    dcc.Graph(id='closeness_centrality', style={'flex': '1'})
                                ])
                            ]),
                            dcc.Tab(label='Shortest Paths', value='shortest-paths', children=[
                                html.Div(style={'display':'flex','alignItems':'center','gap':'10px','padding':'10px'}, children=[
                                    html.Label("Source Residue:"),
                                    dcc.Dropdown(id='source_residue', options=[{'label': r,'value': r} for r in first_res_list],value=None,style={'width':'150px'}),
                                    html.Label("Target Residue:"),
                                    dcc.Dropdown(id='target_residue', options=[{'label': r,'value': r} for r in first_res_list],value=None,style={'width':'150px'}),
                                    html.Button('Find', id='find_paths_btn', n_clicks=0)
                                ]),
                                dash_table.DataTable(
                                    id='paths_table',
                                    columns=[{'name':'Path','id':'Path'},{'name':'Length','id':'Length'}],
                                    data=[],
                                    style_table={'height':'55vh','overflowY':'auto'},
                                    style_cell={'textAlign':'left','whiteSpace':'normal','height':'auto'}
                                )
                            ])
                        ])
                    ])
                ])
            ])
        ]),
        # Right Panel
        html.Div(style={'width': '35%', 'border': '1px solid #CCC', 'padding': '10px','boxSizing':'border-box'}, children=[
            html.H3("3D Molecular Viewer"),
            dash_molstar.MolstarViewer(id='viewer', data=initial_traj, layout={'modelIndex': frame_min}, style={'width': '100%','height':'80%'}),
            html.Div(style={'paddingTop':'10px'}, children=[
                html.Label("Frame:"),
                dcc.Slider(id='frame_slider', min=frame_min, max=frame_max, step=1, value=frame_min,
                           marks={i: str(i) for i in range(frame_min, frame_max+1, max(1,(frame_max-frame_min)//10))},
                           tooltip={'always_visible':True,'placement':'top'})
            ])
        ])
    ])
])

# Callbacks
# Pairwise & Viewer
@app.callback(
    Output('pair_energy_graph','figure'),
    Output('second_residue_table','data'),
    Output('viewer','selection'),
    Output('viewer','focus'),
    Output('viewer','frame'),
    Output('second_residue_table','selected_rows'),
    Input('first_residue_table','selected_rows'),
    Input('second_residue_table','selected_rows'),
    Input('frame_slider','value'),
    State('second_residue_table','data')
)
def update_interface(sel1, sel2, selected_frame, second_data):
    fig = go.Figure(); seldata=None; focusdata=None
    # İlk seçim
    if not sel1:
        return fig, [], no_update, no_update, selected_frame, []
    first = first_res_list[sel1[0]]
    # İkinci tablo
    if not sel2:
        filt = total_df[(total_df['res1']==first)|(total_df['res2']==first)]
        others = [r for r in pd.concat([filt['res1'],filt['res2']]).unique() if r!=first]
        table=[]
        for r in others:
            p1,p2=f"{first}-{r}",f"{r}-{first}"
            vals=total_long[(total_long['Pair']==p1)|(total_long['Pair']==p2)]['Energy']
            ie=round(vals.mean(),3) if not vals.empty else 0
            table.append({'Residue':r,'IE':ie})
        return fig, table, no_update, no_update, selected_frame, []
    # Enerji grafiği ve nokta
    second=second_data[sel2[0]]['Residue']
    p1,p2=f"{first}-{second}",f"{second}-{first}" 
    df_line=total_long[(total_long['Pair']==p1)|(total_long['Pair']==p2)]
    fig.add_trace(go.Scatter(x=df_line['Frame'],y=df_line['Energy'],mode='lines+markers',marker=dict(size=6,opacity=0.5)))
    if selected_frame in df_line['Frame'].astype(int).values:
        e0=df_line[df_line['Frame'].astype(int)==selected_frame]['Energy'].values[0]
        fig.add_trace(go.Scatter(x=[selected_frame],y=[e0],mode='markers',marker=dict(color='red',size=12)))
    fig.update_layout(hovermode='x unified',title=f"Energies for {first}-{second}",xaxis_title='Frame',yaxis_title='Energy')
    # Mol* seçim
    try:
        r1,c1=first.split('_')[0][3:],first.split('_')[1]
        r2,c2=second.split('_')[0][3:],second.split('_')[1]
        t1=molstar_helper.get_targets(c1,r1); t2=molstar_helper.get_targets(c2,r2)
        seldata=molstar_helper.get_selection([t1,t2],select=True,add=False)
        focusdata=molstar_helper.get_focus([t1,t2],analyse=True)
    except:
        seldata,focusdata=no_update,no_update
    return fig, second_data, seldata, focusdata, selected_frame, sel2

# Interaction Matrix
@app.callback(
    Output('matrix_heatmap','figure'),
    Input('frame_slider','value')
)
def update_energy_matrix(frame_value):
    frame_col=str(frame_value)
    if frame_col not in total_df.columns:
        return go.Figure()
    df=total_df[['res1','res2',frame_col]].copy()
    df.columns=['res1','res2','energy']
    df['res1']='E'+df['res1'].str.replace('_','')
    df['res2']='E'+df['res2'].str.replace('_','')
    residues=sorted(set(df['res1']).union(df['res2']))
    matrix_df=pd.DataFrame(0,index=residues,columns=residues)
    for _,row in df.iterrows():
        matrix_df.loc[row['res1'],row['res2']]=row['energy']
        matrix_df.loc[row['res2'],row['res1']]=row['energy']
    fig=go.Figure(data=go.Heatmap(z=matrix_df.values,x=matrix_df.columns,y=matrix_df.index,colorscale='BrBG',zmid=0,zmin=-7,zmax=7,colorbar=dict(title='Energy (kcal/mol)')))
    fig.update_layout(title=f'Interaction Energy Matrix (Frame {frame_value})',xaxis_title='Residue',yaxis_title='Residue',xaxis={'tickangle':45,'automargin':True},yaxis={'automargin':True},margin=dict(l=50,r=50,t=50,b=100),font=dict(size=10))
    return fig

# Network Metrics
@app.callback(
    Output('degree_centrality','figure'),Output('betweenness_centrality','figure'),Output('closeness_centrality','figure'),
    Input('update_network_btn','n_clicks'),Input('frame_slider','value'),Input('include_covalent_edges','value'),Input('energy_cutoff','value')
)
def update_network(n_clicks,frame,include_cov,cutoff):
    G=build_graph(frame,include_cov,cutoff)
    deg=dict(G.degree()); btw=nx.betweenness_centrality(G); clo=nx.closeness_centrality(G)
    def mk(data,title): items=sorted(data.items(),key=lambda x:x[1],reverse=True); nodes,vals=zip(*items); fig=go.Figure(go.Bar(x=vals,y=nodes,orientation='h')); fig.update_layout(title=title,margin=dict(l=150,r=20,t=40,b=40)); return fig
    return mk(deg,'Degree'),mk(btw,'Betweenness centrality'),mk(clo,'Closeness centrality')

# Shortest Paths
@app.callback(
    Output('paths_table','data'),
    Input('find_paths_btn','n_clicks'),Input('frame_slider','value'),Input('include_covalent_edges','value'),Input('energy_cutoff','value'),
    State('source_residue','value'),State('target_residue','value')
)
def find_paths(n_clicks,frame,include_cov,cutoff,source,target):
    if n_clicks<1 or not source or not target or source==target:
        return []
    G=build_graph(frame,include_cov,cutoff)
    try:
        paths_gen=nx.shortest_simple_paths(G,source,target,weight='weight')
        out=[]
        for i,path in enumerate(paths_gen):
            if i>=10: break
            length=sum(G[u][v]['weight'] for u,v in zip(path[:-1],path[1:]))
            out.append({'Path':'-'.join(path),'Length':round(length,6)})
        return out
    except:
        return []

if __name__ == '__main__':
    app.run(debug=True, port=8051)
