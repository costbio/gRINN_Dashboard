from dash import Dash, dcc, html, Input, Output, dash_table
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np


total_df = pd.read_csv("/home/tugba/grinn_workflow_results_new_0/energies_intEnTotal.csv")
elec_df = pd.read_csv("/home/tugba/grinn_workflow_results_new_0/energies_intEnElec.csv")
vdw_df = pd.read_csv("/home/tugba/grinn_workflow_results_new_0/energies_intEnVdW.csv")

total_df = total_df.rename(columns={'Unnamed: 0': 'Pair'})
elec_df = elec_df.rename(columns={'Unnamed: 0': 'Pair'})
vdw_df = vdw_df.rename(columns={'Unnamed: 0': 'Pair'})

# Verileri uzun formatlandı (melt)
total_long = total_df.melt(id_vars=['Pair'], var_name='Frame', value_name='Energy')
elec_long = elec_df.melt(id_vars=['Pair'], var_name='Frame', value_name='Electrostatic Energy')
vdw_long = vdw_df.melt(id_vars=['Pair'], var_name='Frame', value_name='Van der Waals Energy')

pairs_list = total_df['Pair'].unique()    #tekrar kontrol et..
first_amino_acids = sorted(set([p.split('-')[0] for p in pairs_list]))
second_amino_acids = sorted(set([p.split('-')[1] for p in pairs_list]))


# Boş (default) figürler
fig_total = go.Figure()
fig_elec = go.Figure()
fig_vdw = go.Figure()

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
            'padding': 10                     #değiştirilebilir?!!!
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
    
    # Enerji grafik alanı
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


# RENKLENDİRME

def create_colored_line_figure(x_data, y_data, title, yaxis_title="Energy"):
    """
    x_data: Frame veya zaman (numerik olması tercih edilir)
    y_data: İlgili enerji değerleri (Total, Elec, VdW), numerik seri
    title: Grafik başlığı
    yaxis_title: Y ekseni başlığı
    """

    # Her ihtimale karşı sayısal dönüşümler ve NaN temizleme     TEKRAR KONTROL ET!
    x_data = pd.to_numeric(x_data, errors='coerce')
    y_data = pd.to_numeric(y_data, errors='coerce')
    mask = ~(x_data.isna() | y_data.isna())
    x_data = x_data[mask].values
    y_data = y_data[mask].values

    fig = go.Figure()

    if len(x_data) < 2:
        # Veri yok veya tek satır varsa çizgi çizemeyiz
        fig.update_layout(title=title, xaxis_title='Frame', yaxis_title=yaxis_title)
        return fig

    cmin, cmax = float(y_data.min()), float(y_data.max())

    # Her (n-1) segment için ayrı iz ekleyerek renk geçişi sağlıyoruz
    for i in range(len(x_data) - 1):
        # Segmentin X ve Y değerleri
        seg_x = [x_data[i], x_data[i+1]]
        seg_y = [y_data[i], y_data[i+1]]

        # Rengi bu segmentin ortalama y değerine göre belirliyoruz
        avg_val = (y_data[i] + y_data[i+1]) / 2.0
        # 0-1 aralığına normalizasyon
        if cmax != cmin:
            t = (avg_val - cmin) / (cmax - cmin)
        else:
            t = 0.5  # Tüm değerler aynıysa ortalama renk
        t = max(0, min(t, 1))  # emniyet

        # Renk skalasından (Viridis) bir renk seçelim
        color = px.colors.sample_colorscale("Viridis", t)[0]

        fig.add_trace(
            go.Scatter(
                x=seg_x,
                y=seg_y,
                mode='lines',
                line=dict(color=color),
                showlegend=False
            )
        )

    # Renk çubuğu göstermek için DUMMY iz oluşturdum,  veriler geç güncellebilir..
    # Bu iz görünmez olacak, ama colorbar gösterecek
    dummy_y = np.linspace(cmin, cmax, 50)
    dummy_x = np.zeros_like(dummy_y)
    fig.add_trace(
        go.Scatter(
            x=dummy_x,
            y=dummy_y,
            mode='markers',
            marker=dict(
                color=dummy_y,
                colorscale='Viridis',
                cmin=cmin,
                cmax=cmax,
                showscale=True,
                colorbar=dict(title=yaxis_title)
            ),
            opacity=0,
            showlegend=False
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title='Frame',
        yaxis_title=yaxis_title
    )
    return fig



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

#TEKRAR KONTROL ET!
def update_pair_energy_graph(selected_rows_first, selected_rows_second, second_table_data):
    # İlk amino asit seçilmemişse -> tüm 2. amino asitler, boş grafik
    if not selected_rows_first:
        return go.Figure(), go.Figure(), go.Figure(), [{'Pair': aa} for aa in second_amino_acids]
    
    first_amino_acid = first_amino_acids[selected_rows_first[0]]
    
    # 2. tabloyu filtrele
    filtered_second_options = [{'Pair': aa} for aa in second_amino_acids
                               if f"{first_amino_acid}-{aa}" in total_df['Pair'].values]
    
    # 2. amino asit seçilmemişse boş grafik dönsün ama grafik güncellensin
    if not selected_rows_second:
        return go.Figure(), go.Figure(), go.Figure(), filtered_second_options
    
    second_amino_acid = second_table_data[selected_rows_second[0]]['Pair']
    pair = f"{first_amino_acid}-{second_amino_acid}"
    
    if pair not in total_df['Pair'].values:
        return go.Figure(), go.Figure(), go.Figure(), filtered_second_options
    
    # Seçilen pair için veri filtreleme yapılsın
    total_long_2plot = total_long.query(f"Pair == '{pair}'")
    elec_long_2plot = elec_long.query(f"Pair == '{pair}'")
    vdw_long_2plot = vdw_long.query(f"Pair == '{pair}'")

    # Renkli çizgi grafikleri
    fig_total = create_colored_line_figure(
        x_data=total_long_2plot['Frame'],
        y_data=total_long_2plot['Energy'],
        title=f"Total Interaction Energy Time Series for Pair {pair}",
        yaxis_title="Energy"
    )
    fig_elec = create_colored_line_figure(
        x_data=elec_long_2plot['Frame'],
        y_data=elec_long_2plot['Electrostatic Energy'],
        title=f"Electrostatic Energy Time Series for Pair {pair}",
        yaxis_title="Electrostatic Energy"
    )
    fig_vdw = create_colored_line_figure(
        x_data=vdw_long_2plot['Frame'],
        y_data=vdw_long_2plot['Van der Waals Energy'],
        title=f"Van der Waals Energy Time Series for Pair {pair}",
        yaxis_title="Van der Waals Energy"
    )
    
    return fig_total, fig_elec, fig_vdw, filtered_second_options


if __name__ == "__main__":
    app.run_server(debug=True, port="8052")
