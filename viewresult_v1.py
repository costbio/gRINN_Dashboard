import dash, os
from dash import dcc, html, dash_table, Input, Output
import dash_bio as dashbio
import dash_bio.utils.ngl_parser as ngl_parser
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash_bio.utils import PdbParser

# File paths
root      = os.path.dirname(os.path.abspath(__file__))
root      = os.path.join(root, 'test_data', 'prot_lig_1')
pdb_path  = os.path.join(root, 'system_dry.pdb')
total_csv = os.path.join(root, 'energies_intEnTotal.csv')
elec_csv  = os.path.join(root, 'energies_intEnElec.csv')
vdw_csv   = os.path.join(root, 'energies_intEnVdW.csv')

# Read CSVs and prepare energy data
total_df = pd.read_csv(total_csv)
elec_df  = pd.read_csv(elec_csv)
vdw_df   = pd.read_csv(vdw_csv)
for df in (total_df, elec_df, vdw_df):
    df['Pair'] = df['res1'] + '-' + df['res2']

cols2drop = [
    'Unnamed: 0','res1_index','res2_index','res1_chain','res2_chain',
    'res1_resnum','res2_resnum','res1_resname','res2_resname'
]
total_long = total_df.drop(columns=cols2drop+['res1','res2']).melt(
    id_vars=['Pair'], var_name='Frame', value_name='Energy'
)
total_long['Type'] = 'total'
elec_long  = elec_df.drop(columns=cols2drop+['res1','res2']).melt(
    id_vars=['Pair'], var_name='Frame', value_name='Energy'
)
elec_long['Type']  = 'elec'
vdw_long   = vdw_df.drop(columns=cols2drop+['res1','res2']).melt(
    id_vars=['Pair'], var_name='Frame', value_name='Energy'
)
vdw_long['Type']   = 'vdw'
all_energies_long = pd.concat([total_long, elec_long, vdw_long])

first_res_list  = total_df['res1'].unique()
second_res_list = total_df['res2'].unique()

# Parse PDB for chain info
parser     = PdbParser(pdb_path)
model_data = parser.mol3d_data()
atom_df    = pd.DataFrame(model_data['atoms'])

# Build NGL data_list: one entry per chain, colored by chain
chain_ids = atom_df['chain'].unique().tolist()
palette   = px.colors.qualitative.Plotly
base_name = os.path.splitext(os.path.basename(pdb_path))[0]  # "system_dry"

data_list = []
for i, chain in enumerate(chain_ids):
    pdb_id = f"{base_name}.{chain}"  # e.g. "system_dry.A"
    data_list.append(
        ngl_parser.get_data(
            data_path=root + '/',
            pdb_id=pdb_id,
            color=palette[i % len(palette)],
            reset_view=True,
            local=True
        )
    )

# molStyles for NGL viewer: remove axes+box, chain-specific colors
initial_molstyles = {
    "representations": ["cartoon"],  # axes+box removed
    "chosenAtomsColor": "#FFFFFF",
    "chosenAtomsRadius": 1,
    "molSpacingXaxis": 0,
    "sideByside": False
}

# Initialize Dash app
app = dash.Dash(__name__)

# Layout
app.layout = html.Div([
    html.H1("gRINN Workflow Results", style={'textAlign': 'center'}),
    html.Div([

        # Left panel: residue selection & energy plot
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
                                    columns=[{'name': 'Residue', 'id': 'Residue'}],
                                    data=[{'Residue': r} for r in first_res_list],
                                    style_table={'height': '800px', 'overflowY': 'auto', 'width': '120px'},
                                    row_selectable='single',
                                )
                            ], style={'marginLeft': '20px'}),

                            html.Div([
                                html.H4("Select Second Residue"),
                                dash_table.DataTable(
                                    id='second_residue_table',
                                    columns=[{'name': 'Residue', 'id': 'Residue'}],
                                    data=[{'Residue': r} for r in second_res_list],
                                    style_table={'height': '800px', 'overflowY': 'auto', 'width': '120px'},
                                    row_selectable='single',
                                )
                            ], style={'marginRight': '20px'}),

                            html.Div([
                                html.H4("Total Interaction Energy"),
                                dcc.Graph(
                                    id="pair_energy_graph",
                                    figure=go.Figure(),
                                    style={'width': '900px', 'height': '600px'}
                                )
                            ], style={'marginBottom': '20px'})
                        ], style={'display': 'flex', 'flexDirection': 'row', 'marginBottom': '20px'})
                    ]),
                    dcc.Tab(label="Interaction Energy Matrix", value="tab-matrix", children=[html.Div("Matrix...")]),
                    dcc.Tab(label="Interaction Energy Correlations", value="tab-correlations", children=[html.Div("Correlations...")]),
                    dcc.Tab(label="Residue Correlation Matrix", value="tab-residue-matrix", children=[html.Div("Residue Matrix...")]),
                    dcc.Tab(label="Network Analysis", value="tab-network", children=[html.Div("Network...")]),
                ],
                style={'width': '100%', 'height': '100%'},
                colors={'border': '#CCCCCC', 'primary': '#119DFF', 'background': '#F9F9F9'}
            )
        ], style={'width': '55%', 'border': '1px solid #CCC', 'padding': '10px', 'boxSizing': 'border-box'}),

        # Right panel: NGL molecular viewer
        html.Div([
            html.H3("3D Molecular Viewer"),
            dashbio.NglMoleculeViewer(
                id='molecule3d_viewer',
                data=data_list,
                molStyles=initial_molstyles,
                height='800px'
            )
        ], style={'width': '45%', 'border': '1px solid #CCC', 'padding': '10px'})
    ], style={'display': 'flex', 'flexDirection': 'row'})
])

# Callback: update energy plot & second-residue table
@app.callback(
    [Output("pair_energy_graph", "figure"),
     Output("second_residue_table", "data")],
    [Input("first_residue_table", "selected_rows"),
     Input("second_residue_table", "selected_rows"),
     Input("second_residue_table", "data")]
)
def update_results(sel1, sel2, data2):
    if not sel1:
        return go.Figure(), []
    first = first_res_list[sel1[0]]
    filt = total_df[(total_df['res1'] == first) | (total_df['res2'] == first)]
    seconds = pd.concat([filt['res1'], filt['res2']]).unique().tolist()
    seconds = [r for r in seconds if r != first]
    table2 = [{'Residue': r} for r in seconds]
    if not sel2:
        return go.Figure(), table2
    second = data2[sel2[0]]['Residue']
    pair1 = f"{first}-{second}"; pair2 = f"{second}-{first}"
    if pair1 not in total_df['Pair'].values and pair2 not in total_df['Pair'].values:
        return go.Figure(), table2
    df = all_energies_long.query("Pair == @pair1 or Pair == @pair2")
    fig = px.line(df, x="Frame", y="Energy", color="Type", title=f"Energies for {pair1}")
    return fig, table2

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
