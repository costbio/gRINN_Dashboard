import dash, os
from dash import dcc, html, dash_table, Input, Output
import dash_bio as dashbio
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash_bio.utils import PdbParser,  create_mol3d_style

root = os.path.dirname(os.path.abspath(__file__))
root = os.path.join(root,'test_data','prot_lig_1')
pdb_path = os.path.join(root,'system_dry.pdb')
total_csv = os.path.join(root,'energies_intEnTotal.csv')
elec_csv  = os.path.join(root,'energies_intEnElec.csv')
vdw_csv   = os.path.join(root,'energies_intEnVdW.csv')


parser = PdbParser(pdb_path)
model_data = parser.mol3d_data()

total_df = pd.read_csv(total_csv)
elec_df  = pd.read_csv(elec_csv)
vdw_df   = pd.read_csv(vdw_csv)

## TO BE MODIFIED LATER ON
total_df['Pair_indices'] = total_df['res1'] + '-' + total_df['res2']
elec_df['Pair_indices']  = elec_df['res1'] + '-' + elec_df['res2']
vdw_df['Pair_indices']   = vdw_df['res1'] + '-' + vdw_df['res2']

if 'Pair_indices' in total_df.columns:
    total_df = total_df.rename(columns={'Pair_indices': 'Pair'})
if 'Pair_indices' in elec_df.columns:
    elec_df  = elec_df.rename(columns={'Pair_indices': 'Pair'})
if 'Pair_indices' in vdw_df.columns:
    vdw_df   = vdw_df.rename(columns={'Pair_indices': 'Pair'})

#pairs_list = total_df['Pair'].unique()
first_res_list  = total_df['res1'].unique()
second_res_list  = total_df['res2'].unique()

cols2drop = ['Unnamed: 0','res1_index','res2_index','res1_chain','res2_chain','res1_resnum','res2_resnum','res1_resname','res2_resname','res1','res2']

# Drop cols2drop from total_df
total_df_clean = total_df.drop(columns=cols2drop)
elec_df_clean = elec_df.drop(columns=cols2drop)
vdw_df_clean = vdw_df.drop(columns=cols2drop)

total_long = total_df_clean.melt(id_vars=['Pair'], var_name='Frame', value_name='Energy')
total_long['Type'] = 'total'

elec_long = elec_df_clean.melt(id_vars=['Pair'], var_name='Frame', value_name='Energy')
elec_long['Type'] = 'elec'

vdw_long = vdw_df_clean.melt(id_vars=['Pair'], var_name='Frame', value_name='Energy')
vdw_long['Type'] = 'vdw'

all_energies_long = pd.concat([total_long, elec_long, vdw_long])

app = dash.Dash(__name__)

# protein structure
initial_styles = create_mol3d_style(model_data['atoms'], visualization_type= 'cartoon', color_element='residue')

app.layout = html.Div([
    html.H1("gRINN Workflow Results", style={'textAlign': 'center'}),

    html.Div([
        # (%55): Sekmeler + Pairwise Energies
        html.Div([
            dcc.Tabs(
                id="main-tabs",
                value="tab-pairwise",
                children=[
                    dcc.Tab(label="Pairwise Energies", value="tab-pairwise", children=[
                        html.Div([
                            html.Div([
                                html.H4("Select First Residue"),
                                dash_table.DataTable(
                                    id='first_residue_table',
                                    columns=[ 
                                        {'name': 'Residue', 'id': 'Residue'}
                                    ],
                                    data=[{'Residue': res} for res in first_res_list],
                                    style_table={
                                        'height': '800px',
                                        'overflowY': 'auto',
                                        'width': '120px'
                                    },
                                    row_selectable='single',
                                ),
                            ], style={'marginLeft': '20px'}),
                            html.Div([
                                html.H4("Select Second Residue"),
                                dash_table.DataTable(
                                    id='second_residue_table',
                                    columns=[{'name': 'Residue', 'id': 'Residue'}],
                                    data=[{'Residue': res} for res in second_res_list],
                                    style_table={
                                        'height': '800px',
                                        'overflowY': 'auto',
                                        'width': '120px'
                                    },
                                    row_selectable='single',
                                ),
                            ], style={'marginRight': '20px'}),
                            html.Div([
                                html.H4("Total Interaction Energy"),
                                dcc.Graph(
                                    id="pair_energy_graph",
                                    figure=go.Figure(),
                                    style={'width': '700px', 'height': '400px'}
                                )
                            ], style={'marginBottom': '20px'}),
                        ],style={'display': 'flex', 'flexDirection': 'row', 'marginBottom': '20px'}),
                ]),
                    dcc.Tab(label="Interaction Energy Matrix", value="tab-matrix", children=[
                        html.Div("Interaction Energy Matrix içeriği burada...")
                    ]),
                    dcc.Tab(label="Interaction Energy Correlations", value="tab-correlations", children=[
                        html.Div("Interaction Energy Correlations içeriği burada...")
                    ]),
                    dcc.Tab(label="Residue Correlation Matrix", value="tab-residue-matrix", children=[
                        html.Div("Residue Correlation Matrix içeriği burada...")
                    ]),
                    dcc.Tab(label="Network Analysis", value="tab-network", children=[
                        html.Div("Network Analysis içeriği burada...")
                    ]),
                ],
                style={'width': '100%', 'height': '100%'},
                colors={
                    'border': '#CCCCCC',
                    'primary': '#119DFF',
                    'background': '#F9F9F9',
                }
            )
        ],
        style={'width': '55%', 'border': '1px solid #CCC', 'padding': '10px', 'boxSizing': 'border-box'}), 
        # ( %45): 3D Molecular Viewer
        html.Div([
            html.H3("3D Molecular Viewer"),
            dashbio.Molecule3dViewer( 
                id='molecule3d_viewer',
                modelData=model_data,
                styles=initial_styles,  # protein strucutre
                style={'width': '100%', 'height': '800px'},
                selectionType='residue',
                backgroundColor='#FFFFFF',
                zoomTo=True  # protein genişlemiyor????
            ),
        ],
        style={'width': '45%', 'border': '1px solid #CCC', 'padding': '10px', 'marginLeft': '5px', 'boxSizing': 'border-box'}),
    ],
    style={'display': 'flex', 'flexDirection': 'row'})
]) 


@app.callback(
    [
        Output("pair_energy_graph", "figure"),
        #Output("electrostatic_energy_graph", "figure"),
        #Output("vdw_energy_graph", "figure"),
        Output("second_residue_table", "data"),
        Output("molecule3d_viewer", "selectedAtomIds"),
    ],
    [
        Input("first_residue_table", "selected_rows"),
        Input("second_residue_table", "selected_rows"),
        Input("second_residue_table", "data"),
    ]
)
def update_results(selected_rows_first, selected_rows_second, second_table_data):
    if not selected_rows_first:
        return (
            #go.Figure(), go.Figure(), go.Figure(),
            go.Figure(),
            [{'Residue': r} for r in []],
            []
        )

    first_residue = first_res_list[selected_rows_first[0]]

    # Filter total_df to include only rows where first_residue is found in either res1 or res2
    filtered_total = total_df[(total_df['res1'] == first_residue) | (total_df['res2'] == first_residue)]

    # Find all unique res1 and res2 in this filtered dataframe
    second_res_list = pd.concat([filtered_total['res1'], filtered_total['res2']]).unique()

    # First residue should not be in second_res_list
    second_res_list = [res for res in second_res_list if res != first_residue]

    # Filter second_table_data to include only rows where the Residue is in the second_res_list
    filtered_second_options = [{'Residue': res} for res in second_res_list]

    if not selected_rows_second:
        #return (go.Figure(), go.Figure(), go.Figure(), filtered_second_options, [])
        return (go.Figure(), filtered_second_options, [])

    second_residue = second_table_data[selected_rows_second[0]]['Residue']
    pair1 = f"{first_residue}-{second_residue}"
    pair2 = f"{second_residue}-{first_residue}"

    if (pair1 not in total_df['Pair'].values) and (pair2 not in total_df['Pair'].values):
        #return (go.Figure(), go.Figure(), go.Figure(), filtered_second_options, [])
        return (go.Figure(), filtered_second_options, [])

    energies_data = all_energies_long.query("Pair == @pair1 or Pair == @pair2")

    pair_string = f"{first_residue}-{second_residue}"
    fig_total = px.line(energies_data, x="Frame", y="Energy", title=f"Interaction Energies for {pair_string}", color='Type')

    #atom indexleri bul..
    selected_atoms = []
    for i, atom in enumerate(model_data['atoms']):
        if atom['residue_name'] == first_residue or atom['residue_name'] == second_residue:
            selected_atoms.append(i)

    return (fig_total, filtered_second_options, selected_atoms)
    #return (fig_total, fig_elec, fig_vdw, filtered_second_options, selected_atoms)


if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
