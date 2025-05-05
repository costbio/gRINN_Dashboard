import dash, os
from dash import dcc, html, dash_table, Input, Output, State
import dash_bio as dashbio
import dash_bio.utils.ngl_parser as ngl_parser
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash_bio.utils import PdbParser

root      = os.path.dirname(os.path.abspath(__file__))
root      = os.path.join(root, 'test_data', 'prot_lig_1')
pdb_path  = os.path.join(root, 'system_dry.pdb')
total_csv = os.path.join(root, 'energies_intEnTotal.csv')


total_df = pd.read_csv(total_csv)
total_df['Pair'] = total_df['res1'] + '-' + total_df['res2']

cols2drop = [
    'Unnamed: 0','res1_index','res2_index','res1_chain','res2_chain',
    'res1_resnum','res2_resnum','res1_resname','res2_resname'
]

total_long = total_df.drop(columns=cols2drop+['res1','res2']).melt(
    id_vars=['Pair'], var_name='Frame', value_name='Energy'
)
total_long['Type'] = 'total'

first_res_list  = total_df['res1'].unique()
second_res_list = total_df['res2'].unique()


parser     = PdbParser(pdb_path)
model_data = parser.mol3d_data()

# used the ngl package to color according to the protein chains.....?
base_name = os.path.splitext(os.path.basename(pdb_path))[0]
data_list = [
    ngl_parser.get_data(
        data_path  = root + "/",
        pdb_id     = base_name,
        color      = "lightgray",
        reset_view = True,
        local      = True
    )
]

initial_molstyles = {
    "representations":   ["cartoon"],
    "colorScheme":       "chainindex",
    "chosenAtomsColor":  "#FFFFFF",
    "chosenAtomsRadius": 1
}


app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("gRINN Workflow Results", style={'textAlign': 'center'}),
    html.Div([

        
        html.Div([
            dcc.Tabs(
                id="main-tabs",
                value="tab-pairwise",
                style={'width': '100%', 'height': '100%'},
                children=[

                    dcc.Tab(
                        label="Pairwise Energies",
                        value="tab-pairwise",
                        children=html.Div(
                            style={
                                'display': 'flex',
                                'flexDirection': 'row',
                                'height': '100%',
                                'overflow': 'hidden'
                            },
                            children=[

                                # 3 tablo bloğu
                                html.Div(
                                    style={
                                        'display': 'flex',
                                        'flexDirection': 'row',
                                        'flex': '1',
                                        'overflowY': 'auto'
                                    },
                                    children=[
                                        html.Div([
                                            html.H4("Select First Residue"),
                                            dash_table.DataTable(
                                                id='first_residue_table',
                                                columns=[{'name':'Residue','id':'Residue'}],
                                                data=[{'Residue': r} for r in first_res_list],
                                                row_selectable='single',
                                                style_table={
                                                    'height': 'calc(100vh - 300px)',
                                                    'overflowY': 'auto',
                                                    'width': '120px'
                                                }
                                            )
                                        ], style={'marginRight':'20px'}),

                                        html.Div([
                                            html.H4("Select Second Residue"),
                                            dash_table.DataTable(
                                                id='second_residue_table',
                                                columns=[{'name':'Residue','id':'Residue'}],
                                                data=[{'Residue': r} for r in second_res_list],
                                                row_selectable='single',
                                                style_table={
                                                    'height': 'calc(100vh - 300px)',
                                                    'overflowY': 'auto',
                                                    'width': '120px'
                                                }
                                            )
                                        ], style={'marginRight':'20px'}),

                                        html.Div([
                                            html.H4("IE (kcal/mol)"),
                                            dash_table.DataTable(
                                                id='pairwise_summary_table',
                                                columns=[
                                                    {'name':'Residue','id':'Residue'},
                                                    {'name':'IE [kcal/mol]','id':'IE'}
                                                ],
                                                style_table={
                                                    'height': 'calc(100vh - 300px)',
                                                    'overflowY': 'auto',
                                                    'width': '200px'
                                                }
                                            )
                                        ])
                                    ]
                                ),

                            
                                html.Div([
                                    html.H4("Total Interaction Energy"),
                                    dcc.Graph(
                                        id="pair_energy_graph",
                                        figure=go.Figure(),
                                        style={'height':'calc(100vh - 200px)'}
                                    )
                                ], style={
                                    'flex': '2',
                                    'paddingLeft': '20px',
                                    'overflow': 'hidden'
                                })

                            ]
                        )
                    ),

                    
                    dcc.Tab(label="Interaction Energy Matrix", value="tab-matrix", children=[html.Div("Matrix...")]),
                    dcc.Tab(label="Interaction Energy Correlations", value="tab-correlations", children=[html.Div("Correlations...")]),
                    dcc.Tab(label="Residue Correlation Matrix", value="tab-residue-matrix", children=[html.Div("Residue Matrix...")]),
                    dcc.Tab(label="Network Analysis", value="tab-network", children=[html.Div("Network... ")]),
                ]
            )
        ], style={
            'width': '65%',
            'height': '100vh',
            'border': '1px solid #CCC',
            'boxSizing': 'border-box',
            'padding': '10px'
        }),

        
html.Div([
    html.H3("3D Molecular Viewer"),
    html.Div(  
        dashbio.NglMoleculeViewer(
            id='molecule3d_viewer',
            data=data_list,
            molStyles=initial_molstyles,
            width='100%',
            height='100%'
        ),
        style={
            'display': 'flex',
            'justifyContent': 'center',
            'alignItems': 'center',
            'height': 'calc(100vh - 40px)'  # because of header  
        }
    )
], style={
            'width': '35%',
            'height': '100vh',
            'border': '1px solid #CCC',
            'boxSizing': 'border-box',
            'padding': '10px',
            'overflow': 'hidden'
        })

    ], style={'display': 'flex', 'flexDirection': 'row', 'marginTop': '10px'})
])

@app.callback(
    [
        Output("pair_energy_graph", "figure"),
        Output("second_residue_table", "data"),
        Output("pairwise_summary_table", "data"),
    ],
    [
        Input("first_residue_table", "selected_rows"),
        Input("second_residue_table", "selected_rows"),
    ],
    [State("second_residue_table", "data")]
)
def update_results(sel1, sel2, second_data):
    if not sel1:
        return go.Figure(), [], []

    first = first_res_list[sel1[0]]
    filt = total_df[(total_df['res1']==first)|(total_df['res2']==first)]
    seconds = pd.concat([filt['res1'],filt['res2']]).unique().tolist()
    seconds = [r for r in seconds if r!=first]
    second_table = [{'Residue': r} for r in seconds]

    summary = []
    for r in seconds:
        p1, p2 = f"{first}-{r}", f"{r}-{first}"
        subset = total_long.query("Pair==@p1 or Pair==@p2")
        ie_mean = pd.to_numeric(subset['Energy'], errors='coerce').mean() if not subset.empty else None
        summary.append({'Residue': r, 'IE': round(ie_mean,3) if ie_mean is not None else ''})

    if sel2:
        second = second_data[sel2[0]]['Residue']
        p1, p2 = f"{first}-{second}", f"{second}-{first}"
        df_line = total_long.query("Pair==@p1 or Pair==@p2")
        fig = px.line(df_line, x='Frame', y='Energy', title=f"Energies for {p1}")
    else:
        fig = go.Figure()

    return fig, second_table, summary

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
