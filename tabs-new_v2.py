from dash import Dash, dcc, html, Input, Output, dash_table
import pandas as pd
import plotly.express as px

# Veri Yükleme
total_df = pd.read_csv("/home/tugba/grinn_workflow_results_new_0/energies_intEnTotal.csv")
elec_df = pd.read_csv("/home/tugba/grinn_workflow_results_new_0/energies_intEnElec.csv")
vdw_df = pd.read_csv("/home/tugba/grinn_workflow_results_new_0/energies_intEnVdW.csv")

# "Unnamed: 0" sütununu 'Pair' olarak adlandıralım
total_df = total_df.rename(columns={'Unnamed: 0': 'Pair'})
elec_df = elec_df.rename(columns={'Unnamed: 0': 'Pair'})
vdw_df = vdw_df.rename(columns={'Unnamed: 0': 'Pair'})

# Verileri uzun formatlandı (melt)  bunu tekrar arş.
total_long = total_df.melt(id_vars=['Pair'], var_name='Frame', value_name='Energy')
elec_long = elec_df.melt(id_vars=['Pair'], var_name='Frame', value_name='Electrostatic Energy')
vdw_long = vdw_df.melt(id_vars=['Pair'], var_name='Frame', value_name='Van der Waals Energy')

# Tablodaki amino asitleri '-' işareti ile ayırarak liste oluşturuyoruz
pairs_list = total_df['Pair'].unique()
first_amino_acids = sorted(set([pair.split('-')[0] for pair in pairs_list]))
second_amino_acids = sorted(set([pair.split('-')[1] for pair in pairs_list]))

# Boş (default) grafikler
fig_total = px.line(title="Total Interaction Energy Time Series")
fig_elec = px.line(title="Electrostatic Energy Time Series")
fig_vdw = px.line(title="Van der Waals Energy Time Series")

app = Dash()

app.layout = html.Div([
    html.H1("get Residue Interaction Energies and Networks Update Program"),
    
    html.Div([
        # İlk amino asit tablo
        html.Div([
            html.H3("Select First Amino Acid"),
            dash_table.DataTable(
                id='first_amino_acid_table',
                columns=[{'name': 'Pair', 'id': 'Pair'}],
                data=[{'Pair': aa} for aa in first_amino_acids],
                style_table={'height': '300px', 'overflowY': 'auto', 'width': '200px'},
                row_selectable='single',
            ),
        ], style={
            'width': '30%',
            'display': 'inline-block',
            'verticalAlign': 'top',
            'padding': 10
        }),
        
        # İkinci amino asit tablo
        html.Div([
            html.H3("Select Second Amino Acid"),
            dash_table.DataTable(
                id='second_amino_acid_table',
                columns=[{'name': 'Pair', 'id': 'Pair'}],
                data=[{'Pair': aa} for aa in second_amino_acids],
                style_table={'height': '300px', 'overflowY': 'auto', 'width': '200px'},
                row_selectable='single',
            ),
        ], style={
            'width': '30%',
            'display': 'inline-block',
            'verticalAlign': 'top',
            'padding': 10
        }),
    ]),
    
    # Enerji grafiklerinin gösterileceği alan
    html.Div([
        html.Div([ 
            html.H3("Total Interaction Energy"),
            dcc.Graph(id="pair_energy_graph", figure=fig_total)
        ], style={'padding': 10, 'flex': 1, 'width': '60%', 'display': 'inline-block'}),
        
        html.Div([ 
            html.H3("Electrostatic Energy Time Series"),
            dcc.Graph(id="electrostatic_energy_graph", figure=fig_elec)
        ], style={'padding': 10, 'flex': 1, 'width': '60%', 'display': 'inline-block'}),
        
        html.Div([ 
            html.H3("Van der Waals Energy Time Series"),
            dcc.Graph(id="vdw_energy_graph", figure=fig_vdw)
        ], style={'padding': 10, 'flex': 1, 'width': '60%', 'display': 'inline-block'}),
    ], style={'display': 'flex', 'flexDirection': 'row'}),
])


@app.callback(
    [
        Output("pair_energy_graph", "figure"),
        Output("electrostatic_energy_graph", "figure"),
        Output("vdw_energy_graph", "figure"),
        Output("second_amino_acid_table", "data"),
    ],
    [
        Input("first_amino_acid_table", "selected_rows"),
        Input("second_amino_acid_table", "selected_rows"),
        Input("second_amino_acid_table", "data"),
    ]
)
def update_pair_energy_graph(selected_rows_first, selected_rows_second, second_table_data):
    # Eğer ilk amino asit seçilmemişse, tüm ikinci amino asitleri göster, grafikler boş
    if not selected_rows_first:
        return px.line(), px.line(), px.line(), [{'Pair': aa} for aa in second_amino_acids]
    
    first_amino_acid = first_amino_acids[selected_rows_first[0]]
    
    # Filtrelenmiş ikinci amino asitler
    filtered_second_options = [{'Pair': aa} for aa in second_amino_acids
                               if f"{first_amino_acid}-{aa}" in total_df['Pair'].values]
    
    # İkinci amino asit henüz seçilmediyse, grafikler boş kalsın ama tablo güncellensin
    if not selected_rows_second:
        return px.line(), px.line(), px.line(), filtered_second_options
    
    second_amino_acid = second_table_data[selected_rows_second[0]]['Pair']
    pair = f"{first_amino_acid}-{second_amino_acid}"
    
    if pair not in total_df['Pair'].values:
        return px.line(), px.line(), px.line(), filtered_second_options
    
    # Seçilen pair için veriyi filtrele
    total_long_2plot = total_long.query(f"Pair == '{pair}'")
    elec_long_2plot = elec_long.query(f"Pair == '{pair}'")
    vdw_long_2plot = vdw_long.query(f"Pair == '{pair}'")
    
    # 1) Total Interaction Energy grafiği (renk skalası)
    fig_total = px.scatter(
        total_long_2plot,
        x='Frame',
        y='Energy',
        color='Energy',                     # Enerji değeri renk parametresi
        color_continuous_scale='Viridis',   # Renk skalası (isteğe göre değiştirilebilir)
        title=f"Total Interaction Energy Time Series for Pair {pair}"
    )
    fig_total.update_traces(mode='lines+markers')
    
    # 2) Electrostatic Energy grafiği (renk skalası)
    fig_elec = px.scatter(
        elec_long_2plot,
        x='Frame',
        y='Electrostatic Energy',
        color='Electrostatic Energy',
        color_continuous_scale='Viridis',
        title=f"Electrostatic Energy Time Series for Pair {pair}"
    )
    fig_elec.update_traces(mode='lines+markers')
    
    # 3) Van der Waals Energy grafiği (renk skalası)
    fig_vdw = px.scatter(
        vdw_long_2plot,
        x='Frame',
        y='Van der Waals Energy',
        color='Van der Waals Energy',
        color_continuous_scale='Viridis',
        title=f"Van der Waals Energy Time Series for Pair {pair}"
    )
    fig_vdw.update_traces(mode='lines+markers')
    
    return fig_total, fig_elec, fig_vdw, filtered_second_options


if __name__ == "__main__":
    app.run_server(debug=True, port="4051")
