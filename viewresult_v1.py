import os
import dash
from dash import Dash, dcc, html, dash_table, Input, Output, State, callback_context, no_update
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import dash_molstar
from dash_molstar.utils import molstar_helper
from dash_molstar.utils.representations import Representation

# File paths
data_dir = os.path.join(os.path.dirname(__file__), 'test_data', 'prot_lig_1')
pdb_path  = os.path.join(data_dir, 'system_dry.pdb')
total_csv = os.path.join(data_dir, 'energies_intEnTotal.csv')
traj_xtc  = os.path.join(data_dir, 'traj.xtc')

# Load energy data
total_df = pd.read_csv(total_csv)
total_df['Pair'] = total_df['res1'] + '-' + total_df['res2']
cols2drop = [
    'Unnamed: 0','res1_index','res2_index','res1_chain','res2_chain',
    'res1_resnum','res2_resnum','res1_resname','res2_resname'
]
total_long = (
    total_df
    .drop(columns=cols2drop + ['res1','res2'])
    .melt(id_vars=['Pair'], var_name='Frame', value_name='Energy')
)
first_res_list  = total_df['res1'].unique()
second_res_list = total_df['res2'].unique()

# Molstar configuration
cartoon = Representation(type='cartoon', color='uniform')
cartoon.set_color_params({'value': 0xD3D3D3})
chainA = molstar_helper.get_targets(chain='A')
component = molstar_helper.create_component(
    label='Protein', targets=[chainA], representation=cartoon
)

# Prepare trajectory data
topo   = molstar_helper.parse_molecule(pdb_path, component=component)
coords = molstar_helper.parse_coordinate(traj_xtc)
traj_data = molstar_helper.get_trajectory(topo, coords)

# Initialize Dash app
app = Dash(__name__)
app.layout = html.Div([
    html.H1("gRINN Workflow Results", style={'textAlign':'center'}),
    html.Div(
        style={'display':'flex','height':'100vh','gap':'10px'},
        children=[
            # Left panel: Pairwise Energies
            html.Div(
                style={
                    'width':'65%', 'border':'1px solid #CCC', 'padding':'10px',
                    'boxSizing':'border-box'
                },
                children=[
                    dcc.Tabs(
                        id='main-tabs', value='tab-pairwise', children=[
                            dcc.Tab(
                                label='Pairwise Energies', value='tab-pairwise', children=[
                                    html.Div(
                                        style={'display':'flex','height':'calc(100vh - 50px)','gap':'5px'},
                                        children=[
                                            # First residue selection with pagination
                                            html.Div(
                                                style={'flex':'1','overflowY':'auto'},
                                                children=[
                                                    html.H4("Select First Residue"),
                                                    dash_table.DataTable(
                                                        id='first_residue_table',
                                                        columns=[{'name':'Residue','id':'Residue'}],
                                                        data=[{'Residue': r} for r in first_res_list],
                                                        row_selectable='single',
                                                        page_action='native',
                                                        page_current=0,
                                                        page_size=25,
                                                        style_table={
                                                            'height':'calc(100vh - 200px)',
                                                            'width':'120px'
                                                        }
                                                    )
                                                ]
                                            ),
                                            # Second residue + IE with pagination
                                            html.Div(
                                                style={'flex':'1','overflowY':'auto'},
                                                children=[
                                                    html.H4("Select Second Residue & IE"),
                                                    dash_table.DataTable(
                                                        id='second_residue_table',
                                                        columns=[
                                                            {'name':'Residue','id':'Residue'},
                                                            {'name':'IE [kcal/mol]','id':'IE'}
                                                        ],
                                                        data=[],
                                                        row_selectable='single',
                                                        page_action='native',
                                                        page_current=0,
                                                        page_size=15,
                                                        style_table={
                                                            'height':'calc(100vh - 200px)',
                                                            'width':'200px'
                                                        }
                                                    )
                                                ]
                                            ),
                                            # Energy graph
                                            html.Div(
                                                style={'flex':'2','paddingLeft':'20px'},
                                                children=[
                                                    html.H4("Total Interaction Energy"),
                                                    dcc.Graph(
                                                        id='pair_energy_graph',
                                                        style={'height':'calc(100vh - 100px)'}
                                                    )
                                                ]
                                            )
                                        ]
                                    )
                                ]
                            ),
                            dcc.Tab(label='Interaction Energy Matrix', value='tab-matrix', children=[html.Div("Matrix...")]),
                            dcc.Tab(label='Network Analysis', value='tab-network', children=[html.Div("Network...")])
                        ]
                    )
                ]
            ),
            # Right panel: 3D Viewer
            html.Div(
                style={
                    'width':'35%', 'border':'1px solid #CCC', 'padding':'10px',
                    'boxSizing':'border-box'
                },
                children=[
                    html.H3("3D Molecular Viewer"),
                    dash_molstar.MolstarViewer(
                        id='viewer',
                        data=traj_data,
                        style={'width':'100%','height':'90%'}
                    )
                ]
            )
        ]
    )
])

# Unified callback: update table, graph, viewer focus/selection,
# and clear second-residue selection when first table changes.
@app.callback(
    [
        Output('pair_energy_graph','figure'),
        Output('second_residue_table','data'),
        Output('viewer','selection'),
        Output('viewer','focus'),
        Output('second_residue_table','selected_rows')
    ],
    [
        Input('first_residue_table','selected_rows'),
        Input('second_residue_table','selected_rows')
    ],
    [State('second_residue_table','data')]
)
def update_all(sel1, sel2, second_data):
    ctx = callback_context.triggered[0]['prop_id'].split('.')[0]
    fig = go.Figure()
    seldata = None
    focusdata = None

    # 1) First residue changed -> rebuild second table, clear graph & viewer, reset selection
    if ctx == 'first_residue_table':
        if not sel1:
            return fig, [], no_update, no_update, []
        first = first_res_list[sel1[0]]
        filt = total_df[(total_df['res1']==first)|(total_df['res2']==first)]
        seconds = pd.concat([filt['res1'], filt['res2']]).unique().tolist()
        seconds = [r for r in seconds if r!=first]
        table_data = []
        for r in seconds:
            p1, p2 = f"{first}-{r}", f"{r}-{first}"
            subset = total_long.query("Pair==@p1 or Pair==@p2")
            ie = pd.to_numeric(subset['Energy'], errors='coerce').mean() if not subset.empty else 0
            table_data.append({'Residue': r, 'IE': round(ie,3)})
        return fig, table_data, no_update, no_update, []

    # 2) Second residue changed -> update graph & viewer, keep table and selection
    if ctx == 'second_residue_table' and sel1 and sel2:
        first = first_res_list[sel1[0]]
        second = second_data[sel2[0]]['Residue']
        p1, p2 = f"{first}-{second}", f"{second}-{first}"
        df_line = total_long.query("Pair==@p1 or Pair==@p2")
        fig = px.line(df_line, x='Frame', y='Energy', title=f"Energies for {p1}")
        chain1, resid1 = first.split('_')[-1], first.split('_')[0][3:]
        chain2, resid2 = second.split('_')[-1], second.split('_')[0][3:]
        t1 = molstar_helper.get_targets(chain1, resid1)
        t2 = molstar_helper.get_targets(chain2, resid2)
        seldata   = molstar_helper.get_selection([t1,t2], select=True, add=False)
        focusdata = molstar_helper.get_focus([t1,t2], analyse=True)
        return fig, second_data, seldata, focusdata, sel2

    # 3) Fallback: no updates
    return fig, no_update, None, None, no_update

if __name__ == '__main__':
    app.run(debug=True, port=8051)