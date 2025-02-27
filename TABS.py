from dash import Dash, dcc, html, Input, Output, callback, dash_table
import pandas as pd
import plotly.express as px

total_df = pd.read_csv("/mnt/d/experiments/hla_b39/2022_10_14_bioe_leviathan_backup/T1DM_B39_traj_analysis/gmx_traj_data/B3801_1/grinn_output/energies_intEnTotal.csv")

# "Unnamed: 0" sütununu 'Residue' olarak adlandıralım
total_df = total_df.rename(columns={'Unnamed: 0': 'Pair'})

pair_df = pd.DataFrame(total_df['Pair'])

# Verileri uzun formata (melt) dönüştürelim
total_long = total_df.melt(id_vars=['Pair'], var_name='Frame', value_name='Energy')

total_long_2plot = total_long.query('Pair == "345-298"')

fig_total = px.line(total_long_2plot, x='Frame', y='Energy', title="Total Interaction Energy Time Series")

# Starting the Dash app.
app = Dash()
# Uygulamanın layout'u (yapısı)
app.layout = html.Div([

    html.H1("get Residue Interaction Energies and Networks Update Program"),

    # Sayfa arası geçişler için "Tabs"
    dcc.Tabs(id="tabs", value='tab-1', children=[

        dcc.Tab(label="File Upload", value="tab-1"),  # Tab 1: Dosya Yükleme
        dcc.Tab(label="Pairwise Energies", value="tab-2"),
        dcc.Tab(label="Interactions Energy Matrix", value="tab-3"),
        dcc.Tab(label="Interaction Energy Correlations", value="tab-4"),
        dcc.Tab(label="Residue Correlation Matrix", value="tab-5"),
        dcc.Tab(label="Network Analysis", value="tab-6"),

    ]),

    # Sekme içerikleri burada görünecek
    html.Div(id="tabs-content")

])

@callback(Output('tabs-content', 'children'),
              Input('tabs', 'value'))
def render_tab_content(tab):
    if tab == "tab-1":
        return html.Div([
            html.H3("Upload Simulation Data"),

            dcc.Upload(
                id='upload-data',
                children=html.Button('Upload Simulation Data'),
                multiple=True
            ),
            html.Div(id='output-data-upload')
        ])
    elif tab == "tab-2":
        return html.Div([
                    
        html.Div([
            html.H2("Table"),
            dash_table.DataTable(pair_df.to_dict('records'), id="pair_table")
        ], style={'padding': 10, 'flex': 1}),

        html.Div([
            html.H3("Pairwise Energies Tables"),
            dcc.Graph(id="pair_energy_graph",figure=fig_total)
        ], style={'padding': 10, 'flex': 1})], style={'display': 'flex', 'flexDirection': 'row'}

        )

    elif tab == "tab-3":
        return html.Div([
            html.H3("Interactions Energy Matrix Tables")
        ])
    elif tab == "tab-4":
        return html.Div([
            html.H3("Interaction Energy Correlations")
        ])
    elif tab == "tab-5":
        return html.Div([
            html.H3("Residue Correlation Matrix Tables")
        ])
    elif tab == "tab-6":
        return html.Div([
            html.H3("Network Analysis Tables")
        ])

# Sayfa arası geçiş için callback fonksiyonu
@app.callback(
    Output(component_id="pair_energy_graph", component_property="figure"),
    Input(component_id="pair_table", component_property="active_cell")
)
def update_pair_energy_graph(active_cell):
    print(active_cell)
    if active_cell:
        pair = pair_df.loc[active_cell['row'],'Pair']
    else:
        pair = pair_df.loc[0,'Pair']

        total_long_2plot = total_long.query(f'Pair == {pair}')

        fig_total = px.line(total_long_2plot, x='Frame', y='Energy', title="Total Interaction Energy Time Series")
        return fig_total

# Uygulama başlasın
if __name__ == "__main__":
    app.run_server(debug=True, port = "4050")