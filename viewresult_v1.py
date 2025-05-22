import os
import dash
from dash import Dash, dcc, html, dash_table, Input, Output, State, callback_context, no_update
import pandas as pd
import plotly.graph_objects as go
import dash_molstar
from dash_molstar.utils import molstar_helper
from dash_molstar.utils.representations import Representation

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

# Dummy initial pair selection
selected_pair1 = None
selected_pair2 = None

# Extract residue lists
first_res_list = total_df['res1'].unique()

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

# App layout
app = Dash(__name__)
app.layout = html.Div([
    html.H1("gRINN Workflow Results", style={'textAlign': 'center'}),
    html.Div(
        style={'display': 'flex', 'height': '100vh', 'gap': '5px'},
        children=[
            html.Div(
                style={'width': '65%', 'border': '1px solid #CCC', 'padding': '10px', 'boxSizing': 'border-box'},
                children=[
                    dcc.Tabs(id='main-tabs', value='tab-pairwise', children=[
                        dcc.Tab(label='Pairwise Energies', value='tab-pairwise', children=[
                            html.Div(
                                style={'display': 'flex', 'height': 'calc(100vh - 50px)', 'gap': '2px'},
                                children=[
                                    html.Div(
                                        style={'minWidth': '160px', 'maxWidth': '200px', 'overflowY': 'auto'},
                                        children=[
                                            html.H4("Select First Residue"),
                                            dash_table.DataTable(
                                                id='first_residue_table',
                                                columns=[{'name': 'Residue', 'id': 'Residue'}],
                                                data=[{'Residue': r} for r in first_res_list],
                                                row_selectable='single',
                                                style_table={'height': 'calc(100vh - 200px)', 'overflowY': 'scroll'}
                                            )
                                        ]
                                    ),
                                    html.Div(
                                        style={'width': '200px', 'maxWidth': '200px', 'overflowY': 'auto'},
                                        children=[
                                            html.H4("Select Second Residue & IE"),
                                            dash_table.DataTable(
                                                id='second_residue_table',
                                                columns=[
                                                    {'name': 'Residue', 'id': 'Residue'},
                                                    {'name': 'IE [kcal/mol]', 'id': 'IE'}
                                                ],
                                                data=[],
                                                row_selectable='single',
                                                style_table={'height': 'calc(100vh - 200px)', 'overflowY': 'scroll'}
                                            )
                                        ]
                                    ),
                                    html.Div(
                                        style={'flex': '2', 'paddingLeft': '10px'},
                                        children=[
                                            dcc.Graph(
                                                id='pair_energy_graph',
                                                config={'displayModeBar': True},
                                                style={'height': 'calc(100vh - 100px)'}
                                            )
                                        ]
                                    )
                                ]
                            )
                        ]),
                        dcc.Tab(label='Interaction Energy Matrix', value='tab-matrix', children=[html.Div("Matrix...")]),
                        dcc.Tab(label='Network Analysis', value='tab-network', children=[html.Div("Network...")])
                    ])
                ]
            ),
            html.Div(
                style={'width': '35%', 'border': '1px solid #CCC', 'padding': '10px', 'boxSizing': 'border-box'},
                children=[
                    html.H3("3D Molecular Viewer"),
                    dash_molstar.MolstarViewer(
                        id='viewer',
                        data=initial_traj,
                        layout={'modelIndex': frame_min},
                        style={'width': '100%', 'height': '80%'}
                    ),
                    html.Div([
                        html.Label("Frame:", style={'marginBottom': '5px'}),
                        dcc.Slider(
                            id='frame_slider',
                            min=frame_min,
                            max=frame_max,
                            step=1,
                            value=frame_min,
                            marks={i: str(i) for i in range(frame_min, frame_max + 1, max(1, (frame_max - frame_min) // 10))},
                            tooltip={'always_visible': True, 'placement': 'top'}
                        )
                    ], style={'paddingTop': '10px'})
                ]
            )
        ]
    )
])

@app.callback(
    Output('pair_energy_graph', 'figure'),
    Output('second_residue_table', 'data'),
    Output('viewer', 'selection'),
    Output('viewer', 'focus'),
    Output('viewer', 'frame'),
    Output('second_residue_table', 'selected_rows'),
    Input('first_residue_table', 'selected_rows'),
    Input('second_residue_table', 'selected_rows'),
    Input('frame_slider', 'drag_value'),
    State('second_residue_table', 'data')
)
def update_interface(sel1, sel2, selected_frame, second_data):
    
    print(sel1, sel2, selected_frame)
    ctx = callback_context.triggered[0]['prop_id'].split('.')[0]
    fig = go.Figure(); seldata = None; focusdata = None
    
    if ctx == 'first_residue_table':
        if not sel1:
            return fig, [], no_update, no_update, []
        first = first_res_list[sel1[0]]
        filt = total_df[(total_df['res1'] == first) | (total_df['res2'] == first)]
        others = [r for r in pd.concat([filt['res1'], filt['res2']]).unique() if r != first]
        table = []
        for r in others:
            p1 = f"{first}-{r}"; p2 = f"{r}-{first}"
            vals = total_long[(total_long['Pair'] == p1) | (total_long['Pair'] == p2)]['Energy']
            ie = round(vals.mean(), 3) if not vals.empty else 0
            table.append({'Residue': r, 'IE': ie})
        return fig, table, no_update, no_update, selected_frame, []

    if (ctx == 'second_residue_table' and sel1 and sel2) or (ctx == 'frame_slider' and sel1 and sel2):
        first = first_res_list[sel1[0]]
        second = second_data[sel2[0]]['Residue']
        p1 = f"{first}-{second}"; p2 = f"{second}-{first}"
        df_line = total_long[(total_long['Pair'] == p1) | (total_long['Pair'] == p2)]

        fig.add_trace(go.Scatter(
            x=df_line['Frame'],
            y=df_line['Energy'],
            mode='lines+markers',
            name='Energy',
            marker=dict(color='blue', size=6, opacity=0.5),
            line=dict(color='blue')
        ))

        try:
            selected_frame = int(selected_frame)
            if selected_frame in df_line['Frame'].astype(int).values:
                energy_at_frame = df_line[df_line['Frame'].astype(int) == selected_frame]['Energy'].values[0]
                fig.add_trace(go.Scatter(
                    x=[selected_frame], y=[energy_at_frame],
                    mode='markers',
                    marker=dict(color='red', size=12),
                    name='Selected Frame'
                ))
        except Exception as e:
            print("Marker error:", e)


        fig.update_layout(
            clickmode='event+select',
            hovermode='x unified',
            title=f"Energies for {first}-{second}",
            xaxis_title='Frame',
            yaxis_title='Energy'
        )

        c1, r1 = first.split('_')[-1], first.split('_')[0][3:]
        c2, r2 = second.split('_')[-1], second.split('_')[0][3:]
        t1 = molstar_helper.get_targets(c1, r1); t2 = molstar_helper.get_targets(c2, r2)
        seldata = molstar_helper.get_selection([t1, t2], select=True, add=False)
        focusdata = molstar_helper.get_focus([t1, t2], analyse=True)
        return fig, second_data, seldata, focusdata, selected_frame, sel2

    return fig, no_update, None, None, selected_frame, no_update


if __name__ == '__main__':
    app.run(debug=True, port=8051)
