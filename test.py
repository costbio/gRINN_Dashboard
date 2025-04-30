import dash, os
from dash import dcc, html, dash_table, Input, Output
import dash_bio as dashbio
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash_bio.utils import PdbParser, create_mol3d_style