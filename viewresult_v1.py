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

# Residue list - sort by residue number to maintain protein sequence order
def sort_residues_by_sequence(residues):
    """Sort residues by their sequence number extracted from residue names like GLY290_A"""
    def extract_residue_number(res_name):
        try:
            # Extract number from residue name like 'GLY290_A'
            parts = res_name.split('_')
            if len(parts) >= 2:
                # Get the number part from the first part (e.g., '290' from 'GLY290')
                import re
                number = re.findall(r'\d+', parts[0])
                if number:
                    return int(number[0])
            return 0
        except:
            return 0
    
    return sorted(residues, key=extract_residue_number)

first_res_list = sort_residues_by_sequence(total_df['res1'].unique())

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

# St. Patrick's Day theme colors
st_patricks_colors = {
    'primary_green': '#228B22',    # Forest Green
    'light_green': '#90EE90',      # Light Green
    'dark_green': '#006400',       # Dark Green
    'gold': '#FFD700',             # Gold
    'orange': '#FF8C00',           # Dark Orange
    'cream': '#F5F5DC',            # Beige
    'white': '#FFFFFF',            # White
    'emerald': '#50C878',          # Emerald
    'clover': '#3CB371'            # Medium Sea Green
}

# Custom CSS styling
app.index_string = '''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>
            @import url('https://fonts.googleapis.com/css2?family=Irish+Grover&family=Roboto:wght@300;400;500;700&display=swap');
            
            body {
                background: linear-gradient(135deg, #228B22 0%, #90EE90 100%);
                font-family: 'Roboto', sans-serif;
                margin: 0;
                padding: 0;
            }
            
            .main-title {
                font-family: 'Roboto', sans-serif;
                font-weight: 700;
                background: linear-gradient(45deg, #FFD700, #FF8C00);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                background-clip: text;
                text-shadow: 1px 1px 2px rgba(0,0,0,0.2);
                font-size: 2.5rem;
                margin: 20px 0;
                text-align: center;
                position: relative;
                letter-spacing: 1px;
            }
            
            .main-title::before {
                content: "🍀";
                position: absolute;
                left: -60px;
                top: 50%;
                transform: translateY(-50%);
                font-size: 2rem;
            }
            
            .main-title::after {
                content: "🍀";
                position: absolute;
                right: -60px;
                top: 50%;
                transform: translateY(-50%);
                font-size: 2rem;
            }
            
            .panel {
                background: rgba(248,255,248,0.5);
                border: 3px solid #FFD700;
                border-radius: 15px;
                box-shadow: 0 8px 32px rgba(0,0,0,0.2);
                backdrop-filter: blur(10px);
            }
            
            .tab-content {
                background: rgba(250,255,250,0.4);
                border-radius: 10px;
                padding: 15px;
                margin: 10px;
            }
            
            .dash-table-container {
                border-radius: 10px;
                overflow: hidden;
                box-shadow: 0 4px 16px rgba(0,0,0,0.1);
                background: rgba(250,255,250,0.4);
            }
            
            .dash-table-container table {
                background: rgba(250,255,250,0.4) !important;
            }
            
            .dash-table-container .dash-cell {
                background-color: rgba(250,255,250,0.4) !important;
                border-bottom: 1px solid #C0E0C0 !important;
            }
            
            .dash-table-container .dash-header {
                background-color: rgba(50,200,120,0.3) !important;
                color: #2F4F2F !important;
                font-weight: 600 !important;
            }
            
            /* Mol* viewer container */
            .mol-viewer-container {
                background: rgba(240,255,240,0.6);
                border: 2px solid #90EE90;
                border-radius: 10px;
                padding: 10px;
                margin: 10px;
            }
            
            .shamrock-bg {
                background: linear-gradient(45deg, #50C878, #3CB371);
                border-radius: 10px;
                padding: 15px;
                margin: 5px;
                box-shadow: 0 4px 16px rgba(0,0,0,0.1);
            }
            
            /* Button hover effects */
            .dash-button:hover {
                transform: translateY(-2px);
                box-shadow: 0 6px 16px rgba(0,0,0,0.3) !important;
                transition: all 0.3s ease;
            }
            
            /* Dropdown styling */
            .Select-control {
                border: 2px solid #FFD700 !important;
                border-radius: 8px !important;
            }
            
            /* Slider styling */
            .rc-slider-track {
                background-color: #FFD700 !important;
            }
            
            .rc-slider-handle {
                border: 2px solid #FFD700 !important;
                background-color: #228B22 !important;
            }
            
            .rc-slider-dot-active {
                border-color: #FFD700 !important;
            }
            
            /* Tab styling */
            .tab-content {
                animation: fadeIn 0.5s ease-in;
            }
            
            @keyframes fadeIn {
                from { opacity: 0; transform: translateY(20px); }
                to { opacity: 1; transform: translateY(0); }
            }
            
            /* Checkbox styling */
            input[type="checkbox"] {
                transform: scale(1.2);
                margin-right: 8px;
            }
            
            /* Scrollbar styling */
            ::-webkit-scrollbar {
                width: 8px;
            }
            
            ::-webkit-scrollbar-track {
                background: rgba(240,255,240,0.3);
                border-radius: 10px;
            }
            
            ::-webkit-scrollbar-thumb {
                background: #228B22;
                border-radius: 10px;
            }
            
            ::-webkit-scrollbar-thumb:hover {
                background: #006400;
            }
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
'''

app.layout = html.Div([
    html.H1("🍀 gRINN Workflow Results 🍀", 
            className="main-title",
            style={
                'textAlign': 'center',
                'color': st_patricks_colors['gold'],
                'fontFamily': 'Roboto, sans-serif',
                'fontWeight': '700',
                'fontSize': '2.5rem',
                'margin': '20px 0',
                'textShadow': '1px 1px 2px rgba(0,0,0,0.2)',
                'letterSpacing': '1px'
            }),
    html.Div(style={
        'display': 'flex', 
        'height': '100vh', 
        'gap': '15px',
        'padding': '20px',
        'background': f'linear-gradient(135deg, {st_patricks_colors["primary_green"]} 0%, {st_patricks_colors["light_green"]} 100%)'
    }, children=[
        # Left Panel: Tabs
        html.Div(className="panel", style={
            'width': '65%', 
            'padding': '20px', 
            'boxSizing': 'border-box',
            'background': f'rgba(255,255,255,0.95)',
            'border': f'3px solid {st_patricks_colors["gold"]}',
            'borderRadius': '15px',
            'boxShadow': '0 8px 32px rgba(0,0,0,0.2)',
            'backdropFilter': 'blur(10px)'
        }, children=[
            dcc.Tabs(id='main-tabs', value='tab-pairwise', 
                     style={
                         'fontFamily': 'Roboto, sans-serif',
                         'fontWeight': '500'
                     },
                     colors={
                         'border': st_patricks_colors['primary_green'],
                         'primary': st_patricks_colors['gold'],
                         'background': st_patricks_colors['cream']
                     }, children=[
                # Pairwise Energies
                dcc.Tab(label='🔗 Pairwise Energies', value='tab-pairwise', children=[
                    html.Div(className="tab-content", style={
                        'display': 'flex', 
                        'height': 'calc(100vh - 150px)', 
                        'gap': '10px',
                        'background': f'rgba(250,255,250,0.4)',
                        'borderRadius': '10px',
                        'padding': '15px',
                        'margin': '10px'
                    }, children=[
                        html.Div(className="shamrock-bg", style={
                            'minWidth': '180px', 
                            'maxWidth': '220px', 
                            'overflowY': 'auto',
                            'background': f'linear-gradient(45deg, {st_patricks_colors["emerald"]}, {st_patricks_colors["clover"]})',
                            'borderRadius': '10px',
                            'padding': '15px',
                            'boxShadow': '0 4px 16px rgba(0,0,0,0.1)'
                        }, children=[
                            html.H4("🍀 Select First Residue", style={
                                'color': 'white',
                                'fontFamily': 'Roboto, sans-serif',
                                'fontWeight': '500',
                                'textAlign': 'center',
                                'textShadow': '1px 1px 2px rgba(0,0,0,0.5)'
                            }),
                            dash_table.DataTable(
                                id='first_residue_table',
                                columns=[{'name': 'Residue', 'id': 'Residue'}],
                                data=[{'Residue': r} for r in first_res_list],
                                row_selectable='single',
                                style_table={
                                    'height': 'calc(100vh - 250px)', 
                                    'overflowY': 'scroll',
                                    'borderRadius': '8px',
                                    'border': f'2px solid {st_patricks_colors["gold"]}',
                                    'backgroundColor': 'rgba(250,255,250,0.4)'
                                },
                                style_header={
                                    'backgroundColor': st_patricks_colors['gold'],
                                    'color': 'white',
                                    'fontWeight': 'bold',
                                    'textAlign': 'center'
                                },
                                style_cell={
                                    'textAlign': 'center',
                                    'fontFamily': 'Roboto, sans-serif',
                                    'fontSize': '14px',
                                    'backgroundColor': 'rgba(250,255,250,0.4)',
                                    'border': '1px solid #C0E0C0'
                                },
                                style_data_conditional=[
                                    {
                                        'if': {'state': 'selected'},
                                        'backgroundColor': st_patricks_colors['light_green'],
                                        'border': f'2px solid {st_patricks_colors["primary_green"]}'
                                    }
                                ]
                            )
                        ]),
                        html.Div(className="shamrock-bg", style={
                            'width': '220px', 
                            'maxWidth': '220px', 
                            'overflowY': 'auto',
                            'background': f'linear-gradient(45deg, {st_patricks_colors["emerald"]}, {st_patricks_colors["clover"]})',
                            'borderRadius': '10px',
                            'padding': '15px',
                            'boxShadow': '0 4px 16px rgba(0,0,0,0.1)'
                        }, children=[
                            html.H4("🎯 Select Second Residue & IE", style={
                                'color': 'white',
                                'fontFamily': 'Roboto, sans-serif',
                                'fontWeight': '500',
                                'textAlign': 'center',
                                'textShadow': '1px 1px 2px rgba(0,0,0,0.5)'
                            }),
                            dash_table.DataTable(
                                id='second_residue_table',
                                columns=[{'name': 'Residue', 'id': 'Residue'},{'name': 'IE [kcal/mol]', 'id': 'IE'}],
                                data=[],
                                row_selectable='single',
                                style_table={
                                    'height': 'calc(100vh - 250px)', 
                                    'overflowY': 'scroll',
                                    'borderRadius': '8px',
                                    'border': f'2px solid {st_patricks_colors["gold"]}',
                                    'backgroundColor': 'rgba(250,255,250,0.4)'
                                },
                                style_header={
                                    'backgroundColor': st_patricks_colors['gold'],
                                    'color': 'white',
                                    'fontWeight': 'bold',
                                    'textAlign': 'center'
                                },
                                style_cell={
                                    'textAlign': 'center',
                                    'fontFamily': 'Roboto, sans-serif',
                                    'fontSize': '14px',
                                    'backgroundColor': 'rgba(250,255,250,0.4)',
                                    'border': '1px solid #C0E0C0'
                                },
                                style_data_conditional=[
                                    {
                                        'if': {'state': 'selected'},
                                        'backgroundColor': st_patricks_colors['light_green'],
                                        'border': f'2px solid {st_patricks_colors["primary_green"]}'
                                    }
                                ]
                            )
                        ]),
                        html.Div(style={
                            'flex': '2', 
                            'paddingLeft': '15px',
                            'background': 'rgba(250,255,250,0.4)',
                            'borderRadius': '10px',
                            'padding': '15px',
                            'border': f'2px solid {st_patricks_colors["gold"]}'
                        }, children=[
                            dcc.Graph(id='pair_energy_graph', style={'height': 'calc(100vh - 180px)'})
                        ])
                    ])
                ]),
                # Interaction Energy Matrix
                dcc.Tab(label='🔥 Interaction Energy Matrix', value='tab-matrix', children=[
                    html.Div(className="tab-content", style={
                        'background': f'rgba(250,255,250,0.4)',
                        'borderRadius': '10px',
                        'padding': '15px',
                        'margin': '10px',
                        'border': f'2px solid {st_patricks_colors["gold"]}'
                    }, children=[
                        dcc.Graph(id='matrix_heatmap', style={'height': 'calc(100vh - 180px)'})
                    ])
                ]),
                # Network Analysis
                dcc.Tab(label='🕸️ Network Analysis', value='tab-network', children=[
                    html.Div(className="tab-content", style={
                        'background': f'rgba(250,255,250,0.4)',
                        'borderRadius': '10px',
                        'padding': '15px',
                        'margin': '10px',
                        'border': f'2px solid {st_patricks_colors["gold"]}'
                    }, children=[
                        # Controls
                        html.Div(style={
                            'display': 'flex', 
                            'alignItems': 'center', 
                            'paddingBottom': '15px',
                            'background': f'linear-gradient(45deg, {st_patricks_colors["emerald"]}, {st_patricks_colors["clover"]})',
                            'borderRadius': '10px',
                            'padding': '15px',
                            'marginBottom': '15px',
                            'boxShadow': '0 4px 16px rgba(0,0,0,0.1)'
                        }, children=[
                            dcc.Checklist(
                                id='include_covalent_edges',
                                options=[{'label': '🔗 Include covalent bonds as edges', 'value': 'include'}],
                                value=['include'],
                                style={
                                    'marginRight': '20px',
                                    'color': 'white',
                                    'fontFamily': 'Roboto, sans-serif',
                                    'fontWeight': '500'
                                }
                            ),
                            html.Label("⚡ Edge addition energy cutoff (kcal/mol): ", style={
                                'color': 'white',
                                'fontFamily': 'Roboto, sans-serif',
                                'fontWeight': '500'
                            }),
                            dcc.Input(
                                id='energy_cutoff', 
                                type='number', 
                                value=1.0, 
                                step=0.1, 
                                style={
                                    'width': '80px', 
                                    'marginRight': '20px',
                                    'borderRadius': '5px',
                                    'border': f'2px solid {st_patricks_colors["gold"]}',
                                    'padding': '5px'
                                }
                            ),
                            html.Button('🔄 Update Network', 
                                       id='update_network_btn', 
                                       n_clicks=0,
                                       style={
                                           'backgroundColor': st_patricks_colors['gold'],
                                           'color': 'white',
                                           'border': 'none',
                                           'borderRadius': '8px',
                                           'padding': '10px 20px',
                                           'fontWeight': 'bold',
                                           'cursor': 'pointer',
                                           'boxShadow': '0 4px 8px rgba(0,0,0,0.2)'
                                       }),
                            html.Button('💾 Export Network to File...', 
                                       id='export_network_btn', 
                                       n_clicks=0, 
                                       style={
                                           'marginLeft': '10px',
                                           'backgroundColor': st_patricks_colors['orange'],
                                           'color': 'white',
                                           'border': 'none',
                                           'borderRadius': '8px',
                                           'padding': '10px 20px',
                                           'fontWeight': 'bold',
                                           'cursor': 'pointer',
                                           'boxShadow': '0 4px 8px rgba(0,0,0,0.2)'
                                       })
                        ]),
                        # Subtabs
                        dcc.Tabs(id='network-tabs', value='residue-metrics', 
                                colors={
                                    'border': st_patricks_colors['primary_green'],
                                    'primary': st_patricks_colors['gold'],
                                    'background': st_patricks_colors['cream']
                                }, children=[
                            dcc.Tab(label='📊 Residue Metrics', value='residue-metrics', children=[
                                html.Div(style={
                                    'display': 'flex', 
                                    'gap': '20px', 
                                    'height': '55vh', 
                                    'overflowY': 'auto',
                                    'padding': '15px'
                                }, children=[
                                    dcc.Graph(id='degree_centrality', style={'flex': '1'}),
                                    dcc.Graph(id='betweenness_centrality', style={'flex': '1'}),
                                    dcc.Graph(id='closeness_centrality', style={'flex': '1'})
                                ])
                            ]),
                            dcc.Tab(label='🛤️ Shortest Paths', value='shortest-paths', children=[
                                html.Div(style={
                                    'display':'flex',
                                    'alignItems':'center',
                                    'gap':'15px',
                                    'padding':'15px',
                                    'background': f'linear-gradient(45deg, {st_patricks_colors["emerald"]}, {st_patricks_colors["clover"]})',
                                    'borderRadius': '10px',
                                    'marginBottom': '15px'
                                }, children=[
                                    html.Label("🎯 Source Residue:", style={
                                        'color': 'white',
                                        'fontFamily': 'Roboto, sans-serif',
                                        'fontWeight': '500'
                                    }),
                                    dcc.Dropdown(
                                        id='source_residue', 
                                        options=[{'label': r,'value': r} for r in first_res_list],
                                        value=None,
                                        style={'width':'150px'}
                                    ),
                                    html.Label("🏁 Target Residue:", style={
                                        'color': 'white',
                                        'fontFamily': 'Roboto, sans-serif',
                                        'fontWeight': '500'
                                    }),
                                    dcc.Dropdown(
                                        id='target_residue', 
                                        options=[{'label': r,'value': r} for r in first_res_list],
                                        value=None,
                                        style={'width':'150px'}
                                    ),
                                    html.Button('🔍 Find', 
                                               id='find_paths_btn', 
                                               n_clicks=0,
                                               style={
                                                   'backgroundColor': st_patricks_colors['gold'],
                                                   'color': 'white',
                                                   'border': 'none',
                                                   'borderRadius': '8px',
                                                   'padding': '10px 20px',
                                                   'fontWeight': 'bold',
                                                   'cursor': 'pointer',
                                                   'boxShadow': '0 4px 8px rgba(0,0,0,0.2)'
                                               })
                                ]),
                                dash_table.DataTable(
                                    id='paths_table',
                                    columns=[{'name':'Path','id':'Path'},{'name':'Length','id':'Length'}],
                                    data=[],
                                    style_table={
                                        'height':'50vh',
                                        'overflowY':'auto',
                                        'borderRadius': '8px',
                                        'border': f'2px solid {st_patricks_colors["gold"]}'
                                    },
                                    style_header={
                                        'backgroundColor': st_patricks_colors['gold'],
                                        'color': 'white',
                                        'fontWeight': 'bold',
                                        'textAlign': 'center'
                                    },
                                    style_cell={
                                        'textAlign':'left',
                                        'whiteSpace':'normal',
                                        'height':'auto',
                                        'fontFamily': 'Roboto, sans-serif',
                                        'fontSize': '14px'
                                    }
                                )
                            ])
                        ])
                    ])
                ])
            ])
        ]),
        # Right Panel
        html.Div(className="panel", style={
            'width': '35%', 
            'padding': '20px',
            'boxSizing':'border-box',
            'background': f'rgba(255,255,255,0.95)',
            'border': f'3px solid {st_patricks_colors["gold"]}',
            'borderRadius': '15px',
            'boxShadow': '0 8px 32px rgba(0,0,0,0.2)',
            'backdropFilter': 'blur(10px)'
        }, children=[
            html.H3("🧬 3D Molecular Viewer", style={
                'color': st_patricks_colors['primary_green'],
                'fontFamily': 'Roboto, sans-serif',
                'fontWeight': '700',
                'textAlign': 'center',
                'marginBottom': '20px',
                'textShadow': '1px 1px 2px rgba(0,0,0,0.1)'
            }),
            html.Div(style={
                'border': f'3px solid {st_patricks_colors["gold"]}',
                'borderRadius': '10px',
                'overflow': 'hidden',
                'boxShadow': '0 4px 16px rgba(0,0,0,0.1)',
                'marginBottom': '20px',
                'background': 'rgba(250,255,250,0.4)'
            }, children=[
                dash_molstar.MolstarViewer(
                    id='viewer', 
                    data=initial_traj, 
                    layout={'modelIndex': frame_min}, 
                    style={'width': '100%','height':'65vh'}
                )
            ]),
            html.Div(className="shamrock-bg", style={
                'paddingTop':'15px',
                'background': f'linear-gradient(45deg, {st_patricks_colors["emerald"]}, {st_patricks_colors["clover"]})',
                'borderRadius': '10px',
                'padding': '15px',
                'boxShadow': '0 4px 16px rgba(0,0,0,0.1)'
            }, children=[
                html.Label("🎬 Frame:", style={
                    'color': 'white',
                    'fontFamily': 'Roboto, sans-serif',
                    'fontWeight': '500',
                    'fontSize': '16px',
                    'textShadow': '1px 1px 2px rgba(0,0,0,0.5)'
                }),
                dcc.Slider(
                    id='frame_slider', 
                    min=frame_min, 
                    max=frame_max, 
                    step=1, 
                    value=frame_min,
                    marks={i: {
                        'label': str(i),
                        'style': {'color': 'white', 'fontWeight': 'bold'}
                    } for i in range(frame_min, frame_max+1, max(1,(frame_max-frame_min)//10))},
                    tooltip={'always_visible':True,'placement':'top'}
                )
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
    fig.add_trace(go.Scatter(
        x=df_line['Frame'],
        y=df_line['Energy'],
        mode='lines+markers',
        marker=dict(size=6, opacity=0.7, color='#228B22'),
        line=dict(color='#228B22', width=3),
        name='Energy'
    ))
    if selected_frame in df_line['Frame'].astype(int).values:
        e0=df_line[df_line['Frame'].astype(int)==selected_frame]['Energy'].values[0]
        fig.add_trace(go.Scatter(
            x=[selected_frame],
            y=[e0],
            mode='markers',
            marker=dict(color='#FFD700', size=15, symbol='diamond', line=dict(color='#FF8C00', width=2)),
            name='Current Frame'
        ))
    fig.update_layout(
        hovermode='x unified',
        title=f"🍀 Energies for {first}-{second}",
        xaxis_title='Frame',
        yaxis_title='Energy (kcal/mol)',
        plot_bgcolor='rgba(240,255,240,0.3)',
        paper_bgcolor='rgba(255,255,255,0.9)',
        font=dict(family='Roboto, sans-serif', size=12, color='#006400'),
        title_font=dict(family='Roboto, sans-serif', size=16, color='#006400'),
        legend=dict(
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='#FFD700',
            borderwidth=2
        )
    )
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
    
    # Get data for the specific frame
    df=total_df[['res1','res2',frame_col]].copy()
    df.columns=['res1','res2','energy']
    
    # Remove any NaN values
    df = df.dropna()
    
    if df.empty:
        return go.Figure()
    
    # Get all residues and sort them by sequence order
    all_residues = set(df['res1']).union(df['res2'])
    
    # Sort residues by sequence number for proper protein order
    def extract_residue_number(res_name):
        try:
            # Extract number from residue name like 'GLY290_A' 
            import re
            number = re.findall(r'\d+', res_name)
            if number:
                return int(number[0])
            return 0
        except:
            return 0
    
    residues = sorted(all_residues, key=extract_residue_number)
    
    # Create matrix with proper indexing - initialize with NaN to distinguish from zero
    matrix_df=pd.DataFrame(float('nan'), index=residues, columns=residues)
    
    # Fill the matrix with energy values
    for _,row in df.iterrows():
        if row['res1'] in residues and row['res2'] in residues:
            energy_val = float(row['energy'])
            matrix_df.loc[row['res1'],row['res2']] = energy_val
            matrix_df.loc[row['res2'],row['res1']] = energy_val
    
    # Set diagonal to 0 (self-interactions)
    for res in residues:
        matrix_df.loc[res, res] = 0.0
    
    # Fill remaining NaN values with 0
    matrix_df = matrix_df.fillna(0.0)
    
    # Create the heatmap
    fig=go.Figure(data=go.Heatmap(
        z=matrix_df.values,
        x=matrix_df.columns.tolist(),
        y=matrix_df.index.tolist(),
        colorscale='RdYlGn_r',  # St. Patrick's Day inspired colorscale
        zmid=0,
        zmin=-7,
        zmax=7,
        showscale=True,
        colorbar=dict(
            title=dict(
                text='Energy (kcal/mol)',
                font=dict(color='#006400', family='Roboto, sans-serif')
            ),
            tickfont=dict(color='#006400', family='Roboto, sans-serif'),
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='#FFD700',
            borderwidth=2
        )
    ))
    
    fig.update_layout(
        title=f'🔥 Interaction Energy Matrix (Frame {frame_value})',
        xaxis_title='🧬 Residue',
        yaxis_title='🧬 Residue',
        xaxis={'tickangle':45,'automargin':True},
        yaxis={'automargin':True},
        margin=dict(l=80,r=50,t=80,b=100),
        font=dict(size=10, family='Roboto, sans-serif', color='#006400'),
        title_font=dict(size=16, family='Roboto, sans-serif', color='#006400'),
        plot_bgcolor='rgba(240,255,240,0.3)',
        paper_bgcolor='rgba(255,255,255,0.9)',
        height=600
    )
    return fig

# Network Metrics
@app.callback(
    Output('degree_centrality','figure'),Output('betweenness_centrality','figure'),Output('closeness_centrality','figure'),
    Input('update_network_btn','n_clicks'),Input('frame_slider','value'),Input('include_covalent_edges','value'),Input('energy_cutoff','value')
)
def update_network(n_clicks,frame,include_cov,cutoff):
    G=build_graph(frame,include_cov,cutoff)
    deg=dict(G.degree()); btw=nx.betweenness_centrality(G); clo=nx.closeness_centrality(G)
    def mk(data,title): 
        items=sorted(data.items(),key=lambda x:x[1],reverse=True)
        nodes,vals=zip(*items)
        fig=go.Figure(go.Bar(
            x=vals,
            y=nodes,
            orientation='h',
            marker=dict(
                color=vals,
                colorscale='Greens',
                line=dict(color='#FFD700', width=1)
            )
        ))
        fig.update_layout(
            title=f'🍀 {title}',
            margin=dict(l=150,r=20,t=50,b=40),
            font=dict(family='Roboto, sans-serif', size=11, color='#006400'),
            title_font=dict(size=14, family='Roboto, sans-serif', color='#006400'),
            plot_bgcolor='rgba(240,255,240,0.3)',
            paper_bgcolor='rgba(255,255,255,0.9)',
            xaxis=dict(
                title_font=dict(color='#006400'),
                tickfont=dict(color='#006400')
            ),
            yaxis=dict(
                title_font=dict(color='#006400'),
                tickfont=dict(color='#006400')
            )
        )
        return fig
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
