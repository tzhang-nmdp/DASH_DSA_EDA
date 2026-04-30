"""
Main entry point for the Donor Specific Antibody (DSA) Exploratory Data Analysis Dashboard.

This multi-page Dash app uses a fixed sidebar navigation to switch between:
  - Page 1: HSA summary distributions with demographic variables
  - Page 2: 4-way interaction heatmaps of demographic variables on HSA distribution
  - Page 3: HLA matching characteristics comparison (DSA vs non-DSA)
  - Page 4: Transplant outcome / survival analyses (under construction)
"""

import dash
from dash import html, dcc
import dash_bootstrap_components as dbc

# Initialize the multi-page Dash app with Bootstrap SPACELAB theme
app = dash.Dash(__name__, use_pages=True, external_stylesheets=[dbc.themes.SPACELAB])

# --- Sidebar Navigation ---
# Dynamically generates nav links from all registered pages
sidebar = dbc.Nav(
    [
        dbc.NavLink(
            [html.Div(page["name"], className="ms-2")],
            href=page["path"],
            active="exact",
        )
        for page in dash.page_registry.values()
    ],
    vertical=True,
    pills=True,
    className="bg-light",
)

# --- Main Layout ---
# Two-column layout: narrow sidebar (left) + page content (right)
app.layout = dbc.Container([
    # Title row
    dbc.Row([
        dbc.Col(html.Div(
            "Donor Specific Antibody Analyses",
            style={'fontSize': 50, 'textAlign': 'center'}
        ))
    ]),

    html.Hr(),

    # Content row: sidebar + page container
    dbc.Row([
        dbc.Col([sidebar], xs=4, sm=4, md=2, lg=2, xl=2, xxl=2),
        dbc.Col([dash.page_container], xs=8, sm=8, md=10, lg=10, xl=10, xxl=10),
    ])
], fluid=True)


if __name__ == "__main__":
    app.run(debug=False)
