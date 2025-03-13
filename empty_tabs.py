
from dash import Dash, dcc, html, Input, Output, callback

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

# Sayfa arası geçiş için callback fonksiyonu
@app.callback(
    Output("tabs-content", "children"),
    Input("tabs", "value")
)
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
            html.H3("Pairwise Energies Tables")
        ])
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
    
    # Veri yükleme callback fonksiyonu
@app.callback(
    Output('output-data-upload', 'children'),
    Input('upload-data', 'filename')
)
def update_output(uploaded_file):
    if uploaded_file:
        return html.Div(f'Uploaded: {uploaded_file}')
    return 'No file uploaded yet.'


# Uygulama başlasın
if __name__ == "__main__":
    app.run_server(debug=True)