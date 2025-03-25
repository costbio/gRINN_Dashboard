import dash
from dash import dcc, html, dash_table, Input, Output
import dash_bio as dashbio
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash_bio.utils import PdbParser


pdb_path = "/home/tugba/grinn_workflow_results_new_0/system_dry.pdb"
total_csv = "/home/tugba/grinn_workflow_results_new_0/energies_intEnTotal.csv"
elec_csv  = "/home/tugba/grinn_workflow_results_new_0/energies_intEnElec.csv"
vdw_csv   = "/home/tugba/grinn_workflow_results_new_0/energies_intEnVdW.csv"

parser = PdbParser(pdb_path)
model_data = parser.mol3d_data()

total_df = pd.read_csv(total_csv)
elec_df  = pd.read_csv(elec_csv)
vdw_df   = pd.read_csv(vdw_csv)

if 'Pair_indices' in total_df.columns:
    total_df = total_df.rename(columns={'Pair_indices': 'Pair'})
if 'Pair_indices' in elec_df.columns:
    elec_df  = elec_df.rename(columns={'Pair_indices': 'Pair'})
if 'Pair_indices' in vdw_df.columns:
    vdw_df   = vdw_df.rename(columns={'Pair_indices': 'Pair'})

pairs_list = total_df['Pair'].unique()
first_amino_acids  = sorted(set([p.split('-')[0] for p in pairs_list]))
second_amino_acids = sorted(set([p.split('-')[1] for p in pairs_list]))

total_long = total_df.melt(id_vars=['Pair'], var_name='Frame', value_name='Total_Energy')
elec_long  = elec_df.melt(id_vars=['Pair'], var_name='Frame', value_name='Elec_Energy')
vdw_long   = vdw_df.melt(id_vars=['Pair'], var_name='Frame', value_name='VdW_Energy')

app = dash.Dash(__name__)

# protein structure
initial_styles = [
    {
        'visualization_type': 'cartoon',
        'color': 'green'
    }
]

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
                            # ---- TABLOLAR ----
                            html.Div([
                                html.Div([
                                    html.H4("Select First Residue"),
                                    dash_table.DataTable(
                                        id='first_residue_table',
                                        columns=[{'name': 'Residue', 'id': 'Residue'}],
                                        data=[{'Residue': res} for res in first_amino_acids],
                                        style_table={
                                            'height': '350px',
                                            'overflowY': 'auto',
                                            'width': '120px'
                                        },
                                        row_selectable='single',
                                    ),
                                ], style={'marginRight': '20px'}),
                                html.Div([
                                    html.H4("Select Second Residue"),
                                    dash_table.DataTable(
                                        id='second_residue_table',
                                        columns=[{'name': 'Residue', 'id': 'Residue'}],
                                        data=[{'Residue': res} for res in second_amino_acids],
                                        style_table={
                                            'height': '350px',
                                            'overflowY': 'auto',
                                            'width': '120px'
                                        },
                                        row_selectable='single',
                                    ),
                                ]),
                            ],
                            style={'display': 'flex', 'flexDirection': 'row', 'marginBottom': '20px'}),

                            # ---- GRAFİKLER ----
                            html.Div([
                                html.Div([
                                    html.H4("Total Interaction Energy"),
                                    dcc.Graph(
                                        id="pair_energy_graph",
                                        figure=go.Figure(),
                                        style={'width': '700px', 'height': '400px'}
                                    )
                                ], style={'marginBottom': '20px'}),
                                html.Div([
                                    html.H4("Electrostatic Energy"),
                                    dcc.Graph(
                                        id="electrostatic_energy_graph",
                                        figure=go.Figure(),
                                        style={'width': '700px', 'height': '400px'}
                                    )
                                ], style={'marginBottom': '20px'}),
                                html.Div([
                                    html.H4("Van der Waals Energy"),
                                    dcc.Graph(
                                        id="vdw_energy_graph",
                                        figure=go.Figure(),
                                        style={'width': '700px', 'height': '400px'}
                                    )
                                ]),
                            ],
                            style={'display': 'flex', 'flexDirection': 'column'}),
                        ],
                        style={'display': 'flex', 'flexDirection': 'column', 'margin': '10px'}),
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
                style={'width': '100%'},
                colors={
                    'border': '#CCCCCC',
                    'primary': '#119DFF',
                    'background': '#F9F9F9'
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
],
style={'fontFamily': 'Arial, sans-serif'})


@app.callback(
    [
        Output("pair_energy_graph", "figure"),
        Output("electrostatic_energy_graph", "figure"),
        Output("vdw_energy_graph", "figure"),
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
            go.Figure(), go.Figure(), go.Figure(),
            [{'Residue': r} for r in second_amino_acids],
            []
        )

    first_residue = first_amino_acids[selected_rows_first[0]]
    filtered_second = [
        r for r in second_amino_acids if f"{first_residue}-{r}" in total_df['Pair'].values
    ]
    filtered_second_options = [{'Residue': res} for res in filtered_second]

    if not selected_rows_second:
        return (go.Figure(), go.Figure(), go.Figure(), filtered_second_options, [])

    second_residue = second_table_data[selected_rows_second[0]]['Residue']
    pair = f"{first_residue}-{second_residue}"

    if pair not in total_df['Pair'].values:
        return (go.Figure(), go.Figure(), go.Figure(), filtered_second_options, [])

    total_data = total_long.query("Pair == @pair")
    elec_data  = elec_long.query("Pair == @pair")
    vdw_data   = vdw_long.query("Pair == @pair")

    fig_total = px.line(total_data, x="Frame", y="Total_Energy", title=f"Total Energy for {pair}")
    fig_elec  = px.line(elec_data,  x="Frame", y="Elec_Energy",  title=f"Electrostatic Energy for {pair}")
    fig_vdw   = px.line(vdw_data,   x="Frame", y="VdW_Energy",   title=f"Van der Waals Energy for {pair}")

    #atom indexleri bul..
    selected_atoms = []
    for i, atom in enumerate(model_data['atoms']):
        if atom['residue_name'] == first_residue or atom['residue_name'] == second_residue:
            selected_atoms.append(i)

    return (fig_total, fig_elec, fig_vdw, filtered_second_options, selected_atoms)


if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
