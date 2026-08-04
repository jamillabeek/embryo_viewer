import dash
from dash import html, dcc, Output, Input, State
from dash.dependencies import ALL, MATCH
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import plotly.express as px
import scanpy as sc
import numpy as np
import json, os, re
import pandas as pd 
from flask_caching import Cache
from dash import dash_table
from dash import Patch


CLUSTER_MARKERS_FILE = "data/cluster_marker_list.csv"
SUBCLUSTER_MARKERS_FILE = "data/subcluster_marker_list.csv"

df_cluster_markers = pd.read_csv(CLUSTER_MARKERS_FILE)
df_subcluster_markers = pd.read_csv(SUBCLUSTER_MARKERS_FILE)


############## define data ###############
H5AD_FILE = "data/adata.h5ad"

CLUSTER_COL = "whole_leiden"
SUBCLUSTER_COL = "sub_leiden"
SECTION_COL = "section"

X_COL = "x_centroid"
Y_COL = "z_centroid"
Z_COL = "y_centroid"

COLORS_FILE = "data/colors.txt"
COLORS_FILE_PREFIX = ""
COLORS_FILE_CLUSTER_COLUMN = "cell cluster"
COLORS_FILE_COLOR_COLUMN = "color code"

MESH_FOLDER = "data/scaffolds"

PANEL_SPLIT = [{"label": "Embryo", "prefix": "e"},{"label": "Placenta", "prefix": "p"},]

# disable split:
#PANEL_SPLIT = None



# ---------- load data ----------
mtx = sc.read_h5ad(H5AD_FILE)


section_labels = sorted(mtx.obs[SECTION_COL].unique())
gene_list = mtx.var_names.tolist()

mtx.obs["_cl"] = mtx.obs[CLUSTER_COL].astype(str)
mtx.obs["_subcl"] = mtx.obs[SUBCLUSTER_COL].astype(str)
mtx.obs["_section"] = mtx.obs[SECTION_COL].astype(str)


def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]

cluster_color_map = {}

if COLORS_FILE:
    cl_col = pd.read_csv(COLORS_FILE, sep="\t")

    cl_col["cluster"] = (
        COLORS_FILE_PREFIX
        + cl_col[COLORS_FILE_CLUSTER_COLUMN].astype(str)
    )

    cluster_color_map = (
        cl_col
        .set_index("cluster")[COLORS_FILE_COLOR_COLUMN]
        .to_dict()
    )
else:
    import colorsys
    n = len(mtx.obs["_cl"].unique())
    cluster_color_map = {
        c: '#{:02x}{:02x}{:02x}'.format(
            *[int(x * 255) for x in colorsys.hsv_to_rgb(i / n, 0.7, 0.9)]
        )
        for i, c in enumerate(sorted(mtx.obs["_cl"].unique()))
    }
cluster_labels = sorted(mtx.obs["_cl"].unique(), key=natural_key)


cluster_colors = {
    c: cluster_color_map.get(c, "#808080")
    for c in cluster_labels
}

cluster_groups = {}

if PANEL_SPLIT:
    for panel in PANEL_SPLIT:
        cluster_groups[panel["label"]] = sorted(
            [c for c in cluster_labels if c.startswith(panel["prefix"])],
            key=natural_key
        )
else:
    cluster_groups["Clusters"] = cluster_labels

cluster_to_group = {}

for group, clusters in cluster_groups.items():
    for c in clusters:
        cluster_to_group[c] = group
    
parent_to_sub = {}

for parent in cluster_labels:
    subs = (
        mtx.obs.loc[
            mtx.obs["_cl"] == parent,
            "_subcl"
        ]
        .replace("", np.nan)
        .dropna()
        .unique()
    )

    if len(subs):
        parent_to_sub[parent] = sorted(map(str, subs), key=natural_key)

parent_map = {sub: parent for parent, subs in parent_to_sub.items() for sub in subs}
subcluster_colors = {sub: cluster_colors.get(parent_map[sub], '#808080') for sub in parent_map}

def build_cluster_panel(label, clusters):
    print(label, len(clusters))
    id = {"type": "toggle-all", "group": label}

    return html.Div([
        html.Div([
            html.Div(label, style={'fontWeight':'bold'}),
            html.Div("Surface",
                     id={'type': 'surface-header', 'group': label},
                     style={'textAlign':'center',
                            'fontWeight':'bold'}),
        ], style={
            'display':'grid',
            'gridTemplateColumns':'auto 34px 18px',
            'columnGap':'2px',
            'marginBottom':'4px'
        }),
        dcc.Checklist(
            id=id,
            options=[{'label': 'All', 'value': 'all'}],
            value=[],
            style={'marginBottom':'6px'}
        ), 

        html.Div(
            control_rows(
                clusters,
                parent_to_sub,
                cluster_colors,
                subcluster_colors,
                mesh_data
            )
        )
    ], style={
        'flex':'1',
        'minWidth':'100px'
    })
# ---------- include scaffolds ----------

mesh_data = {}

def _cluster_key_from_fname(fname: str):
    m = re.match(r"(.+)_mesh\.json$", fname)  # "e_0_mesh.json"
    if m: return m.group(1)
    m = re.match(r"mesh_(.+)\.json$", fname)  # "mesh_e_0.json"
    if m: return m.group(1)
    return None

if os.path.isdir(MESH_FOLDER):
    for fname in os.listdir(MESH_FOLDER):
        if not fname.endswith(".json"): continue
        key = _cluster_key_from_fname(fname)
        if key is None: continue
        try:
            with open(os.path.join(MESH_FOLDER, fname), "r") as f:
                mesh_data[key] = json.load(f)
        except Exception as e:
            print(f"⚠️ Failed to load {fname}: {e}")
print("Loaded meshes for:", sorted(mesh_data.keys()))

# ---------- UI helpers ----------
def mesh_checkbox(id_dict):
    return dcc.Checklist(
        id=id_dict,
        options=[{'label': '', 'value': 'mesh'}],
        value=[],
        style={'transform': 'scale(0.8)', 'margin': '0'}
    )

def control_rows(clusters, parent_to_sub, cluster_colors, subcluster_colors, mesh_data):
    """Build tidy 3-col grid rows: [label+check] [color] [mesh] for parents + subclusters."""
    blocks = []
    for parent in clusters:
        # parent row
        prow = html.Div([
            dcc.Checklist(
                id={'type':'cluster-check','index':parent},
                options=[{'label': parent, 'value': parent}],
                value=[]
            ),
            dcc.Input(
                id={'type':'cluster-color-input','index':parent},
                type='color',
                value=cluster_colors.get(parent,'#808080'),
                style={'width':'40px','height':'24px','padding':'0','border':'none'}
            ),
            mesh_checkbox({'type':'mesh-toggle','index':parent}) if parent in mesh_data else html.Div()
        ], style={'display': 'flex', 'alignItems': 'center', 'gap': '4px'})

        # subcluster rows
        srows = []
        for sub in parent_to_sub.get(parent, []):
            srows.append(
                html.Div([
                    dcc.Checklist(
                        id={'type':'subcluster-check','index':sub},
                        options=[{'label': sub, 'value': sub}],
                        value=[]
                    ),
                    dcc.Input(
                        id={'type':'subcluster-color-input','index':sub},
                        type='color',
                        value=subcluster_colors.get(sub, cluster_colors.get(parent,'#808080')),
                        style={'width':'36px','height':'22px','padding':'0','border':'none'}
                    ),
                    mesh_checkbox({'type':'mesh-toggle','index':sub}) if sub in mesh_data else html.Div()
                ], style={'display': 'flex', 'alignItems': 'center', 'gap': '4px'})
            )
        blocks.append(html.Div([
            prow,
            html.Details([
                html.Summary("subclusters", style={'fontSize': '11px', 'color': '#666', 'cursor': 'pointer', 'paddingLeft': '20px'}),
                html.Div(srows, style={'paddingLeft': '20px'})
            ]) if srows else html.Div()
        ]))    
    return blocks


LOW_OPACITY_KEYS = {"e2", "e5"}
LOW_OPACITY_VALUE = 0.2
FULL_OPACITY_VALUE = 1.0

def build_welcome_figure():
    fig = go.Figure()
    items = sorted(mesh_data.items(), key=lambda kv: kv[0] in LOW_OPACITY_KEYS)
    for key, mesh in items:
        verts = np.array(mesh['verts'])
        faces = np.array(mesh['faces'], dtype=int)
        if faces.size == 0 or verts.size == 0:
            continue
        i, j, k = faces.T
        color = cluster_color_map.get(key, '#808080')
        opacity = LOW_OPACITY_VALUE if key in LOW_OPACITY_KEYS else FULL_OPACITY_VALUE
        fig.add_trace(go.Mesh3d(
            x=verts[:,0], y=verts[:,1], z=verts[:,2],
            i=i, j=j, k=k,
            color=color,
            opacity=opacity,
            flatshading=False,
            showlegend=False, hoverinfo="name",
        ))

    angle = np.radians(20)
    base_x, base_y, base_z = 1.25, 1.25, 1.25
    eye_x = base_x * np.cos(angle) - base_y * np.sin(angle)
    eye_y = base_x * np.sin(angle) + base_y * np.cos(angle)

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            bgcolor='white',
            aspectmode='data',
            camera=dict(eye=dict(x=eye_x, y=eye_y, z=base_z))
        ),
        paper_bgcolor='white',
        margin=dict(l=0, r=0, t=0, b=0),
        uirevision='welcome'
    )
    return fig

# ---------- Dash app layout ----------
app = dash.Dash(__name__)

cache = Cache(app.server, config={'CACHE_TYPE': 'SimpleCache'})


app.layout = html.Div([

    html.H1("3D Cluster Viewer", style={'textAlign': 'center'}),

    html.Div([
        # cluster selection
        html.Div(
            id='sidebar',
            children=[
                html.Div([
                    html.Label("Select up to 3 Genes (for RGB coloring):", style={'fontWeight': 'bold'}),
                    dcc.Dropdown(
                        id='gene_selector',
                        options=[{'label': g, 'value': g} for g in gene_list],
                        multi=True, value=[],
                        placeholder="Pick genes to color cells",
                        style={'width': '100%'}
                    ),
                    html.Label("Gamma:", style={'fontSize': '12px'}),
                    dcc.Slider(id='gamma_slider', min=0.1, max=2.0, step=0.1, value=1.0,
                            tooltip={"placement": "bottom", "always_visible": False}),
                    html.Div(id='tab1-only-sliders', children=[
                        html.Label("Mesh Opacity:", style={'fontSize': '12px'}),
                        dcc.Slider(id='mesh_opacity', min=0.05, max=1.0, step=0.05, value=0.25,
                                tooltip={"placement": "bottom", "always_visible": False}),
                        html.Label("Section spacing (Z Zoom):", style={'fontSize': '12px'}),
                        dcc.Slider(id='z_zoom_slider', min=0.25, max=2.0, step=0.05, value=1.0,
                                tooltip={"placement": "bottom", "always_visible": False}),
                    ]),
                ], style={'paddingBottom': '12px', 'borderBottom': '1px solid #ccc', 'marginBottom': '12px'}),

                html.Div(
                    [build_cluster_panel(label, clusters)
                    for label, clusters in cluster_groups.items()],
                    style={'display': 'flex', 'flexDirection': 'row', 'gap': '8px'}
                ),
            ],
            style={
                'flex': '1', 'display': 'flex', 'flexDirection': 'column',
                'maxHeight': '90vh', 'overflowY': 'auto',
                'paddingLeft': '4px', 'minWidth': '180px'
            }
        ),

        html.Div(
            dcc.Tabs(id='tabs', value='tab-0', children=[
                dcc.Tab(label='Welcome', value='tab-0', children=[
                    html.Div([
                        html.Div([
                            html.H2("3D atlas of Carnegie stage 10 embryo"),
                            html.P("Welcome to our interactive companion viewer for 'A single-cell spatial atlas of Carnegie stage 10 human embryogenesis', [Journal, Year]."),
                            html.P("This viewer allows for exploration of our cell clusters, subclusters and their marker genes in 3D or 2D space. You can select clusters and subclusters to display, adjust their colors, and choose up to 3 genes to visualize using the panel on the right, and you can navigate between the 3D view tab and the sections view (2D) tab. Some clusters in the 3D viewer also have surface meshes that can be toggled on or off, to approach the structure of the intact tissue. The cluster numbers correspond to supplementary table 1 in the publication."),
                            html.P("Please note that the viewer includes a large amount of data, when viewing the full dataset it may take up to a few minutes to load, please wait for the viewer to finish loading before interacting, and be a bit patient when switching tabs or making significant changes."),
                            html.P("For more information or any issues with the viewer, please refer to the paper and the supplementary materials or contact us via xxx@xxx.nl."),
                            ], style={
                            'flex': '1', 'padding': '40px', 'display': 'flex',
                            'flexDirection': 'column', 'justifyContent': 'center'
                        }),

                        html.Div([
                            dcc.Graph(id='welcome_plot', style={'height': '80vh', 'width': '100%'})
                        ], style={'flex': '1', 'display': 'flex', 'alignItems': 'center', 'justifyContent': 'center'})
                    ], style={'display': 'flex', 'flexDirection': 'row', 'minHeight': '80vh'})
                ]),
                #3d viewer tab
                dcc.Tab(label='3D View', value='tab-1', children=[
                    html.Div([
                        html.Div([
                            html.Label("Select Sections to Display:"),
                            dcc.Checklist(
                                id='section_selector',
                                options=[{'label': str(sec), 'value': sec} for sec in section_labels],
                                value=section_labels,
                                inline=True,
                                style={'marginBottom': '10px'}
                            ),
                            html.Label("Marker Size:", style={'fontSize': '12px'}),
                            dcc.Slider(id='size_slider', min=0.1, max=5, step=0.1, value=0.3,
                                tooltip={"placement": "bottom", "always_visible": False}),
                        ], style={'flex': '1', 'minWidth': '280px', 'padding': '0 20px'})
                    ], style={'display': 'flex', 'gap': '20px', 'marginBottom': '16px', 'alignItems': 'flex-start'}),

                    html.Div([
                        html.Div([
                            dcc.Graph(id='cluster_3d_plot', style={'height': '600px', 'width': '100%'})
                        ], style={'flex': '1', 'minWidth': '0', 'marginRight': '8px'}),

                        html.Div([
                            dcc.Graph(id='gene_expression_plot', style={'height': '600px', 'width': '100%'}),
                            html.Div(id='rgb_legend', style={
                                'position': 'absolute', 'top': '40px', 'right': '10px',
                                'zIndex': '10', 'backgroundColor': 'rgba(255,255,255,0.7)',
                                'borderRadius': '8px', 'padding': '6px'
                            })
                        ], style={'position': 'relative', 'flex': '1', 'minWidth': '0', 'marginLeft': '8px'})

                    ], style={'display': 'flex', 'flexDirection': 'row', 'flex': '3', 'minWidth': '0'}),
                ]),
                # sections view tab
                dcc.Tab(label='Sections View', value='tab-2', children=[
                    html.Div([
                        html.Div([
                            html.Label("Select Sections:"),
                            dcc.Checklist(
                                id='section_selector_t2',
                                options=[{'label': str(s), 'value': s} for s in section_labels],
                                value=section_labels,
                                inline=True
                            ),
                            html.Label("Marker Size:", style={'fontSize': '12px'}),
                            dcc.Slider(id='size_slider_t2', min=0.1, max=5, step=0.1, value=1.5,
                                tooltip={"placement": "bottom", "always_visible": False}),
                        ], style={'marginBottom': '12px'}),

                        html.Div(id='section_panels', style={
                            'display': 'grid',
                            'gridTemplateColumns': 'repeat(5, 1fr)',
                            'gap': '8px',
                        }),
                    ])
                ]),

            ]),
            style={'flex': '4', 'minWidth': '0'}
        ),

    ], style={'display': 'flex', 'gap': '8px', 'alignItems': 'flex-start'}),
    # marker tables
    html.Div(id='marker_tables', children=[
        html.Div([
            dcc.Dropdown(
                id='marker_search_t1',
                options=[{'label': g, 'value': g} for g in sorted(df_cluster_markers['names'].unique())],
                multi=True,
                placeholder="Search marker genes..."
            ),
            dash_table.DataTable(
                id='cluster_marker_table',
                columns=[{'name': c, 'id': c} for c in df_cluster_markers.columns],
                data=df_cluster_markers.to_dict('records'),
                page_size=15,
                sort_action='native',
                style_table={'overflowX': 'auto'}
            )
        ], style={'flex': '1'}),

        html.Div([
            dcc.Dropdown(
                id='marker_search_t2',
                options=[{'label': g, 'value': g} for g in sorted(df_subcluster_markers['names'].unique())],
                multi=True,
                placeholder="Search marker genes..."
            ),
            dash_table.DataTable(
                id='subcluster_marker_table',
                columns=[{'name': c, 'id': c} for c in df_subcluster_markers.columns],
                data=df_subcluster_markers.to_dict('records'),
                page_size=15,
                sort_action='native',
                style_table={'overflowX': 'auto'}
            )
        ], style={'flex': '1'}),
    ], style={'display': 'flex', 'gap': '16px', 'marginTop': '20px'}),

    dcc.Store(id='cluster_selector', data=[]),
    dcc.Store(id='cluster_colors_store', data=cluster_colors),
    dcc.Store(id='subcluster_colors_store', data=subcluster_colors),
    dcc.Store(id='camera-store'),
    dcc.Store(id='debounce-store'),
    dcc.Interval(id='debounce-interval', interval=3000, n_intervals=0, disabled=True),

])
## --------- Debugging --------


@app.callback(
    Output('cluster_colors_store', 'data', allow_duplicate=True),
    Input('cluster_colors_store', 'data'),
    prevent_initial_call=True
)
def debug_cluster_colors(data):
    print("cluster_colors_store value:", data)
    return dash.no_update
# @app.callback(
#     Output('cluster_selector', 'data'),
#     Input({'type': 'cluster-check', 'index': ALL}, 'value'),
# )
# def sync_selected_clusters(values):
#     selected_clusters = [val for v in (values or []) for val in (v or [])]
#     return selected_clusters

@app.callback(
    Output('gene_expression_plot', 'figure', allow_duplicate=True),
    Input('gene_expression_plot', 'figure'),
    prevent_initial_call=True
)
def debug_gene_expression(fig):
    print("gene_expression_plot figure updated")
    return dash.no_update
# -----------------------------------------------

#@app.callback(
#    Output('camera-store', 'data'),
#    Input('cluster_3d_plot', 'relayoutData'),
#    Input('gene_expression_plot', 'relayoutData'),
#    prevent_initial_call=True
#)
#def save_camera_state(relayout_left, relayout_right):
#    trigger = getattr(ctx, "triggered_id", None)
#    if not trigger:
#        trig_raw = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None
#        trigger = json.loads(trig_raw) if trig_raw and trig_raw.startswith('{') else trig_raw

#    if trigger == 'cluster_3d_plot' and relayout_left and 'scene.camera' in relayout_left:
#        return relayout_left['scene.camera']

#    if trigger == 'gene_expression_plot' and relayout_right and 'scene.camera' in relayout_right:
#        return relayout_right['scene.camera']

#    raise dash.exceptions.PreventUpdate

@app.callback(
    Output('welcome_plot', 'figure'),
    Input('tabs', 'value'),
)
def render_welcome(active_tab):
    if active_tab != 'tab-0':
        raise PreventUpdate
    return build_welcome_figure()

# write checkbox states to a store and enable interval for debouncing
@app.callback(
    Output('debounce-store', 'data'),
    Output('debounce-interval', 'disabled'),
    Input({'type': 'cluster-check', 'index': ALL}, 'value'),
    prevent_initial_call=True
)
def buffer_cluster_changes(values):
    return values, False

@app.callback(
    Output('cluster_selector', 'data'),
    Output('debounce-interval', 'disabled'),
    Input('debounce-interval', 'n_intervals'),
    State('debounce-store', 'data'),
    prevent_initial_call=True
)
def apply_debounced(_, buffered):
    selected = [val for v in (buffered or []) for val in (v or [])]
    return selected, True

# to hide the surface checkboxes in tab2:
@app.callback(
    Output('sidebar', 'style'),
    Output('tab1-only-sliders', 'style'),
    Output({'type': 'mesh-toggle', 'index': ALL}, 'style'),
    Output({'type': 'surface-header', 'group': ALL}, 'style'),
    Input('tabs', 'value'),
    State({'type': 'mesh-toggle', 'index': ALL}, 'id'),
    State({'type': 'surface-header', 'group': ALL}, 'id'),
)
def toggle_tab_visibility(active_tab, mesh_ids, header_ids):
    hidden = {'display': 'none'}
    visible = {}

    sidebar_style = {
        'flex': '1', 'display': 'flex', 'flexDirection': 'column',
        'maxHeight': '90vh', 'overflowY': 'auto',
        'paddingLeft': '4px', 'minWidth': '180px',
    }
    if active_tab == 'tab-0':
        sidebar_style = {**sidebar_style, 'display': 'none'}

    sliders_style = hidden if active_tab == 'tab-2' else visible
    mesh_styles = [hidden if active_tab == 'tab-2' else visible] * len(mesh_ids)
    header_styles = [hidden if active_tab == 'tab-2' else visible] * len(header_ids)

    return sidebar_style, sliders_style, mesh_styles, header_styles

@app.callback(
    Output('sidebar', 'style'),
    Output('tab1-only-sliders', 'style'),
    Output({'type': 'mesh-toggle', 'index': ALL}, 'style'),
    Output({'type': 'surface-header', 'group': ALL}, 'style'),
    Output('marker_tables', 'style'),
    Input('tabs', 'value'),
    State({'type': 'mesh-toggle', 'index': ALL}, 'id'),
    State({'type': 'surface-header', 'group': ALL}, 'id'),
)
def toggle_tab_visibility(active_tab, mesh_ids, header_ids):
    hidden = {'display': 'none'}
    visible = {}

    sidebar_style = {
        'flex': '1', 'display': 'flex', 'flexDirection': 'column',
        'maxHeight': '90vh', 'overflowY': 'auto',
        'paddingLeft': '4px', 'minWidth': '180px',
    }
    tables_style = {'display': 'flex', 'gap': '16px', 'marginTop': '20px'}

    if active_tab == 'tab-0':
        sidebar_style = {**sidebar_style, 'display': 'none'}
        tables_style = {**tables_style, 'display': 'none'}

    sliders_style = hidden if active_tab == 'tab-2' else visible
    mesh_styles = [hidden if active_tab == 'tab-2' else visible] * len(mesh_ids)
    header_styles = [hidden if active_tab == 'tab-2' else visible] * len(header_ids)

    return sidebar_style, sliders_style, mesh_styles, header_styles, tables_style

# --- Toggle all clusters ---
@app.callback(
    Output({'type': 'cluster-check', 'index': ALL}, 'value'),

    Input({'type': 'toggle-all', 'group': ALL}, 'value'),

    State({'type': 'toggle-all', 'group': ALL}, 'id'),
    State({'type': 'cluster-check', 'index': ALL}, 'id'),
    State({'type': 'cluster-check', 'index': ALL}, 'value'),
    prevent_initial_call=True
)
def toggle_all(group_values, group_ids, cluster_ids, current_vals):

    trigger = ctx.triggered_id

    if not trigger:
        raise PreventUpdate

    if trigger.get("type") != "toggle-all":
        raise PreventUpdate

    target_group = trigger["group"]

    toggle_state = None

    for value, gid in zip(group_values, group_ids):
        if gid["group"] == target_group:
            toggle_state = value
            break

    turn_on = "all" in (toggle_state or [])

    out = list(current_vals or ([[]] * len(cluster_ids)))

    for i, cid in enumerate(cluster_ids):

        cluster = cid["index"]

        if cluster_to_group.get(cluster) == target_group:
            out[i] = [cluster] if turn_on else []

    return out



@app.callback(
    Output({'type': 'subcluster-check', 'index': ALL}, 'value'),
    Input({'type': 'cluster-check', 'index': ALL}, 'value'),
    State({'type': 'cluster-check', 'index': ALL}, 'id'),
    State({'type': 'subcluster-check', 'index': ALL}, 'id'),
    prevent_initial_call=True
)
def parent_toggles_children(parent_values, parent_ids, sub_ids):

    checked_parents = {pid['index']
                       for val, pid in zip(parent_values or [], parent_ids or [])
                       if (val and len(val) > 0)}

    out = []
    for sid in sub_ids or []:
        sub = sid['index']
        parent = parent_map.get(sub)
        if parent in checked_parents:
            out.append([sub])          # selected
        else:
            out.append([])             # deselected
    return out

# ----------------------------------------

# --- Update cluster colors ---
@app.callback(
    Output('cluster_colors_store', 'data'),
    Input({'type': 'cluster-color-input', 'index': ALL}, 'value'),
    State({'type': 'cluster-color-input', 'index': ALL}, 'id'),
    State('cluster_colors_store', 'data')
)
def update_cluster_colors(colors, ids, current_colors):
    base = dict(current_colors or {})
    if colors is None or ids is None:
        return base
    for color, cid in zip(colors, ids):
        if color and cid and 'index' in cid:
            base[cid['index']] = color
    return base

@app.callback(
    Output('subcluster_colors_store', 'data'),
    Input({'type': 'subcluster-color-input', 'index': ALL}, 'value'),
    State({'type': 'subcluster-color-input', 'index': ALL}, 'id'),
    State('subcluster_colors_store', 'data'),
    prevent_initial_call=True
)
def update_subcluster_colors(colors, ids, current):
    base = dict(current or {})
    for color, cid in zip(colors or [], ids or []):
        if color and 'index' in cid:
            base[cid['index']] = color
    return base

try:
    from dash import ctx
except Exception:
    ctx = dash.callback_context

@app.callback(
    Output({'type': 'subcluster-color-input', 'index': ALL}, 'value'),
    Input({'type': 'cluster-color-input', 'index': ALL}, 'value'),
    State({'type': 'cluster-color-input', 'index': ALL}, 'id'),
    State({'type': 'subcluster-color-input', 'index': ALL}, 'id'),
)
def apply_parent_color_to_subs(parent_colors, parent_ids, sub_ids):
    if not parent_colors or not parent_ids:
        return [dash.no_update] * len(sub_ids or [])

    triggered = getattr(ctx, "triggered_id", None)
    if not triggered:
        t = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None
        triggered = json.loads(t) if t else None

    if not triggered or triggered.get('type') != 'cluster-color-input':
        return [dash.no_update] * len(sub_ids or [])

    changed_parent = triggered['index']
    parent_color_map = {pid['index']: col for col, pid in zip(parent_colors, parent_ids)}
    new_color = parent_color_map.get(changed_parent)

    if not new_color:
        return [dash.no_update] * len(sub_ids or [])

    out = []
    for sid in sub_ids or []:
        sub = sid['index']
        if parent_map.get(sub) == changed_parent:
            out.append(new_color)
        else:
            out.append(dash.no_update)
    return out

def rgb_array_to_hex(rgb):
    rgb_uint8 = (np.clip(rgb, 0, 1) * 255).astype(np.uint8)
    return ['#%02x%02x%02x' % (r, g, b) for r, g, b in rgb_uint8]

@cache.memoize(timeout=0)
def _gene_vec(adata, gene):
    """Return 1-D numpy array of expression for 'gene' across all cells."""
    Xg = adata[:, gene].X
    if hasattr(Xg, "toarray"):
        arr = Xg.toarray()
    else:
        arr = np.asarray(Xg)
    return np.ravel(arr)

# filtering the tables
@app.callback(
    Output('cluster_marker_table', 'data'),
    Input('cluster_selector', 'data'),
    Input('marker_search_t1', 'value'),
)
def update_cluster_marker_table(selected_clusters, search_genes):
    df = df_cluster_markers[df_cluster_markers['whole_leiden'].isin(selected_clusters or [])]
    if search_genes:
        df = df[df['names'].isin(search_genes)]
    return df.to_dict('records')

@app.callback(
    Output('subcluster_marker_table', 'data'),
    Input('cluster_selector', 'data'),
    Input('marker_search_t2', 'value'),
)
def update_subcluster_marker_table(selected_clusters, search_genes):
    df = df_subcluster_markers[df_subcluster_markers['whole_leiden'].isin(selected_clusters or [])]
    if search_genes:
        df = df[df['names'].isin(search_genes)]
    return df.to_dict('records')

# lighter color patch

@app.callback(
    Output('cluster_3d_plot', 'figure', allow_duplicate=True),
    Output('gene_expression_plot', 'figure', allow_duplicate=True),
    Input({'type': 'cluster-color-input', 'index': ALL}, 'value'),
    Input({'type': 'subcluster-color-input', 'index': ALL}, 'value'),
    State({'type': 'cluster-color-input', 'index': ALL}, 'id'),
    State({'type': 'subcluster-color-input', 'index': ALL}, 'id'),
    State('cluster_3d_plot', 'figure'),
    State('gene_expression_plot', 'figure'),
    prevent_initial_call=True
)
def patch_colors(parent_colors, sub_colors, parent_ids, sub_ids, fig1_state, fig2_state):
    color_map = {}
    for col, cid in zip(parent_colors or [], parent_ids or []):
        if cid:
            color_map[cid['index']] = col
    for col, cid in zip(sub_colors or [], sub_ids or []):
        if cid:
            color_map[cid['index']] = col

    patch1, patch2 = Patch(), Patch()
    changed = False

    for i, tr in enumerate(fig1_state['data']):
        name = tr.get('name', '')
        for key, color in color_map.items():
            if name.startswith(f'{key} '):
                patch1['data'][i]['marker']['color'] = color
                changed = True

    for i, tr in enumerate(fig2_state['data']):
        name = tr.get('name', '')
        for key, color in color_map.items():
            if name.startswith(f'{key} ') or f' in {key}' in name:
                patch2['data'][i]['marker']['color'] = color
                changed = True

    if not changed:
        raise PreventUpdate
    return patch1, patch2
# ----------------------------------------



# @app.callback(
#     Output('cluster_3d_plot', 'figure'),
#     Output('gene_expression_plot', 'figure'),
#     Output('rgb_legend', 'children'),
#     Input('cluster_selector', 'data'),
#     Input('section_selector', 'value'),
#     Input('gene_selector', 'value'),
#     Input('gamma_slider', 'value'),
#     Input('size_slider', 'value'),
#     Input('mesh_opacity', 'value'),
#     Input('z_zoom_slider', 'value'),
#     Input({'type': 'mesh-toggle', 'index': ALL}, 'value'),
#     Input({'type': 'subcluster-check', 'index': ALL}, 'value'),
#     # Input({'type': 'cluster-color-input', 'index': ALL}, 'value'),
#     # Input({'type': 'subcluster-color-input', 'index': ALL}, 'value'),
#     State('cluster_colors_store', 'data'),
#     State('subcluster_colors_store', 'data'),
#     State({'type': 'mesh-toggle', 'index': ALL}, 'id'),
#     State({'type': 'subcluster-check', 'index': ALL}, 'id'),
#     State({'type': 'cluster-color-input', 'index': ALL}, 'id'),
#     State({'type': 'subcluster-color-input', 'index': ALL}, 'id'),
#     State('camera-store', 'data'),
#     Input('tabs', 'value'),
# )
# def update_figures(
#     selected_clusters, selected_sections, selected_genes,
#     gamma, marker_size, mesh_opacity, z_zoom,
#     mesh_values, sub_values, parent_colors, sub_colors,
#     mesh_ids, sub_ids, parent_color_ids, sub_color_ids, camera_state, active_tab
# ):
@app.callback(
    Output('cluster_3d_plot', 'figure'),
    Output('gene_expression_plot', 'figure'),
    Output('rgb_legend', 'children'),
    Input('cluster_selector', 'data'),
    Input('section_selector', 'value'),
    Input('gene_selector', 'value'),
    Input('gamma_slider', 'value'),
    Input('size_slider', 'value'),
    Input('mesh_opacity', 'value'),
    Input('z_zoom_slider', 'value'),
    Input({'type': 'mesh-toggle', 'index': ALL}, 'value'),
    Input({'type': 'subcluster-check', 'index': ALL}, 'value'),
    State('cluster_colors_store', 'data'),
    State('subcluster_colors_store', 'data'),
    State({'type': 'mesh-toggle', 'index': ALL}, 'id'),
    State({'type': 'subcluster-check', 'index': ALL}, 'id'),
    State('camera-store', 'data'),
    Input('tabs', 'value'),
)
def update_figures(
    selected_clusters, selected_sections, selected_genes,
    gamma, marker_size, mesh_opacity, z_zoom,
    mesh_values, sub_values, parent_colors, sub_colors,
    mesh_ids, sub_ids, camera_state, active_tab
):
    if active_tab != 'tab-1':
        raise PreventUpdate
    triggered = {t["prop_id"] for t in ctx.triggered}
    gene_only_change = all(
        "gene_selector" in t or "gamma_slider" in t
        for t in triggered
    )
    selected_clusters = selected_clusters or []

    mesh_toggle_dict = {
        cid['index']: val
        for cid, val in zip(mesh_ids or [], mesh_values or [])
    }

    # live_cluster_colors = dict(cluster_colors)
    # if parent_colors and parent_color_ids:
    #     for col, cid in zip(parent_colors, parent_color_ids):
    #         if col and cid and 'index' in cid:
    #             live_cluster_colors[cid['index']] = col

    # live_subcluster_colors = dict(subcluster_colors)
    # if sub_colors and sub_color_ids:
    #     for col, cid in zip(sub_colors, sub_color_ids):
    #         if col and cid and 'index' in cid:
    #             live_subcluster_colors[cid['index']] = col
    live_cluster_colors = dict(parent_colors or cluster_colors)
    live_subcluster_colors = dict(sub_colors or subcluster_colors)

    # Which subclusters are currently checked?
    selected_subs = {
        sid['index']
        for val, sid in zip(sub_values or [], sub_ids or [])
        if (val and len(val) > 0)
    }

    if selected_genes is None:
        selected_genes = []
    elif isinstance(selected_genes, str):
        selected_genes = [selected_genes]
    else:
        # It may be a tuple/np array; coerce and drop falsy entries
        selected_genes = [g for g in list(selected_genes) if g]

    # y_factor = z_zoom or 1.0
    # base_mask = mtx.obs[SECTION_COL].isin(selected_sections or section_labels)
    # if selected_clusters:
    #     base_mask &= mtx.obs['_cl'].isin(selected_clusters)

    # y0 = float(np.nanmedian(mtx.obs.loc[base_mask, Y_COL])) if base_mask.any() else 0.0
    # def scale_y(yvals):
    #     return (yvals - y0) * y_factor + y0
    
    y_factor = z_zoom or 1.0
    base_mask = mtx.obs[SECTION_COL].isin(selected_sections or section_labels)
    if selected_clusters:
        base_mask &= mtx.obs['_cl'].isin(selected_clusters)
    y0 = float(np.nanmedian(mtx.obs.loc[base_mask, Y_COL])) if base_mask.any() else 0.0

    def scale_y(arr):
        arr = np.asarray(arr)
        return (arr - y0) * y_factor + y0
    
    # ---- Figure 1: clusters ----
    fig1 = go.Figure()
    for parent in selected_clusters:
        parent_mask = (
            (mtx.obs['_cl'] == parent) &
            (mtx.obs[SECTION_COL].isin(selected_sections))
        )
        all_child_subs = parent_to_sub.get(parent, [])
        child_subs = [s for s in parent_to_sub.get(parent, []) if s in selected_subs]

        if child_subs:
            for sub in child_subs:
                mask = parent_mask & (mtx.obs['_subcl'] == sub)
                fig1.add_trace(go.Scatter3d(
                    x=mtx.obs.loc[mask, X_COL],
                    y=scale_y(mtx.obs.loc[mask, Y_COL]),
                    z=mtx.obs.loc[mask, Z_COL],
                    mode='markers',
                    marker=dict(
                        size=marker_size,
                        color=live_subcluster_colors.get(sub, live_cluster_colors.get(parent, '#808080'))
                    ),
                    name=f'{sub} cells',
                    legendgroup=parent, hoverinfo="name"
                ))

                # optional sub mesh
                if ('mesh' in mesh_toggle_dict.get(sub, [])) and (sub in mesh_data):
                    mesh = mesh_data[sub]
                    verts = np.array(mesh['verts']); faces = np.array(mesh['faces'], dtype=int)
                    if faces.size > 0 and verts.size > 0:
                        verts[:,1] = scale_y(verts[:,1])
                        i, j, k = faces.T
                        fig1.add_trace(go.Mesh3d(
                            x=verts[:,0], y=verts[:,1], z=verts[:,2],
                            i=i, j=j, k=k,
                            color=live_subcluster_colors.get(sub, live_cluster_colors.get(parent, '#808080')),
                            opacity=mesh_opacity,
                            flatshading=False,
                            name=f"{sub} mesh",
                            legendgroup=parent,
                            showlegend=True, hoverinfo="name"
                        ))
        elif not all_child_subs:
            fig1.add_trace(go.Scatter3d(
            x=mtx.obs.loc[parent_mask, X_COL],
            y=scale_y(mtx.obs.loc[parent_mask, Y_COL]),
            z=mtx.obs.loc[parent_mask, Z_COL],
            mode='markers',
            marker=dict(size=marker_size, color=live_cluster_colors.get(parent, '#808080')),
            name=f'{parent} cells',
            legendgroup=parent
        ))
            

        if ('mesh' in mesh_toggle_dict.get(parent, [])) and (parent in mesh_data):
            mesh = mesh_data[parent]
            verts = np.array(mesh['verts']); faces = np.array(mesh['faces'], dtype=int)
            if faces.size > 0 and verts.size > 0:
                verts[:, 1] = scale_y(verts[:, 1])
                i, j, k = faces.T
                fig1.add_trace(go.Mesh3d(
                    x=verts[:,0], y=verts[:,1], z=verts[:,2],
                    i=i, j=j, k=k,
                    color=live_cluster_colors.get(parent, '#808080'),
                    opacity=mesh_opacity,
                    flatshading=False,
                    name=f"{parent} mesh",
                    legendgroup=parent,
                    showlegend=True, hoverinfo="name"
                ))

    print(camera_state)
    fig1.update_layout(
        scene=dict(
            xaxis=dict(title='', showbackground=False, showticklabels=False),
            yaxis=dict(title='', showbackground=False, showticklabels=False),
            zaxis=dict(title='', showbackground=False, showticklabels=False),
            aspectmode='data'
        ),
        margin=dict(l=0,r=0,b=0,t=30),
        showlegend=False,
        uirevision='constant'
    )
    #if camera_state:
    #    fig1.update_layout(scene_camera=camera_state)

    # ---- Figure 2: gene expression ----
    fig2 = go.Figure()
    rgb_legend = []

    # centered Y scaling
    # y_factor = z_zoom or 1.0
    # base_mask = mtx.obs[SECTION_COL].isin(selected_sections or section_labels)
    # if selected_clusters:
    #     base_mask &= mtx.obs['_cl'].isin(selected_clusters)
    # y0 = float(np.nanmedian(mtx.obs.loc[base_mask, Y_COL])) if base_mask.any() else 0.0

    # def scale_y(arr):
    #     arr = np.asarray(arr)
    #     return (arr - y0) * y_factor + y0

    for cluster in selected_clusters:
        mask_cells = (
            (mtx.obs['_cl'] == cluster) &
            (mtx.obs[SECTION_COL].isin(selected_sections))
        )
        if not mask_cells.any():
            continue
        all_child_subs = parent_to_sub.get(cluster, [])
        child_subs = [s for s in all_child_subs if s in selected_subs]

        # Pre-computed arrays (ALL consistently scaled or unscaled)
        x_arr = mtx.obs.loc[mask_cells, X_COL].to_numpy()
        y_arr = scale_y(mtx.obs.loc[mask_cells, Y_COL].to_numpy())   # scaled here
        z_arr = mtx.obs.loc[mask_cells, Z_COL].to_numpy()

        # --- Parent mesh (optional) ---
        if ('mesh' in mesh_toggle_dict.get(cluster, [])) and (cluster in mesh_data):
            mesh = mesh_data[cluster]
            verts = np.array(mesh['verts']).copy()
            faces = np.array(mesh['faces'], dtype=int)
            if faces.size > 0 and verts.size > 0:
                verts[:, 1] = scale_y(verts[:, 1])
                i, j, k = faces.T
                fig2.add_trace(go.Mesh3d(
                    x=verts[:,0], y=verts[:,1], z=verts[:,2],
                    i=i, j=j, k=k,
                    color='lightgrey',
                    opacity=mesh_opacity,
                    flatshading=False,
                    showlegend=False, hoverinfo="name"
                ))

        # --- Subcluster meshes ---
        for sub in child_subs:
            if ('mesh' in mesh_toggle_dict.get(sub, [])) and (sub in mesh_data):
                smesh = mesh_data[sub]
                sverts = np.array(smesh['verts']).copy()
                sfaces = np.array(smesh['faces'], dtype=int)
                if sfaces.size > 0 and sverts.size > 0:
                    sverts[:, 1] = scale_y(sverts[:, 1])
                    i, j, k = sfaces.T
                    fig2.add_trace(go.Mesh3d(
                        x=sverts[:,0], y=sverts[:,1], z=sverts[:,2],
                        i=i, j=j, k=k,
                        color='lightgrey',
                        opacity=mesh_opacity,
                        flatshading=False,
                        name=f"{sub} mesh",
                        showlegend=False, hoverinfo="name"
                    ))

        # 0 genes: color by sub if any selected, else by parent color
        if len(selected_genes) == 0:
            if child_subs:
                for sub in child_subs:
                    sub_mask = mask_cells & (mtx.obs['_subcl'] == sub)
                    fig2.add_trace(go.Scatter3d(
                        x=mtx.obs.loc[sub_mask, X_COL],
                        y=scale_y(mtx.obs.loc[sub_mask, Y_COL]),
                        z=mtx.obs.loc[sub_mask, Z_COL],
                        mode='markers',
                        marker=dict(
                            size=marker_size,
                            color=live_subcluster_colors.get(sub, live_cluster_colors.get(cluster, '#808080')),
                            opacity=0.85
                        ),
                        name=f'{sub} cells', hoverinfo="name"
                    ))
            # else:
            #     fig2.add_trace(go.Scatter3d(
            #         x=x_arr, y=y_arr, z=z_arr,
            #         mode='markers',
            #         marker=dict(
            #             size=marker_size,
            #             color=live_cluster_colors.get(cluster, '#808080'),
            #             opacity=0.85
            #         ),
            #         name=f'{cluster} cells', hoverinfo="name"
            #     ))
            # continue
            elif not all_child_subs:
                fig2.add_trace(go.Scatter3d(
                    x=x_arr, y=y_arr, z=z_arr,
                    mode='markers',
                    marker=dict(size=marker_size, color=live_cluster_colors.get(cluster, '#808080'), opacity=0.85),
                    name=f'{cluster} cells', hoverinfo="name"
                ))
            continue

        # 1 gene: scalar colors
        if len(selected_genes) == 1:
            gene = selected_genes[0]
            if gene in mtx.var_names:
                gv_full = _gene_vec(mtx, gene)
                if child_subs:
                    for sub in child_subs:
                        sub_mask = mask_cells & (mtx.obs['_subcl'] == sub)
                        sub_np = sub_mask.to_numpy()
                        gv = gv_full[sub_np]
                        if gv.size and np.nanmax(gv) > 0:
                            gv = np.log1p(gv); gmin, gmax = float(np.nanmin(gv)), float(np.nanmax(gv))
                            norm = ((gv - gmin) / (gmax - gmin + 1e-9)) ** gamma
                        else:
                            norm = np.zeros_like(gv, dtype=float)
                        fig2.add_trace(go.Scatter3d(
                            x=mtx.obs.loc[sub_mask, X_COL],
                            y=scale_y(mtx.obs.loc[sub_mask, Y_COL]),
                            z=mtx.obs.loc[sub_mask, Z_COL],
                            mode='markers',
                            marker=dict(size=marker_size, color=norm,
                                        colorscale='Viridis', cmin=0, cmax=1, opacity=0.9,
                                        colorbar=dict(title=f"{gene} expr", len=0.5)),
                            name=f"{gene} in {sub}", hoverinfo="name"
                        ))
                elif not all_child_subs:
                    gv = gv_full[mask_cells.to_numpy()]
                    if gv.size and np.nanmax(gv) > 0:
                        gv = np.log1p(gv); gmin, gmax = float(np.nanmin(gv)), float(np.nanmax(gv))
                        norm = ((gv - gmin) / (gmax - gmin + 1e-9)) ** gamma
                    else:
                        norm = np.zeros_like(gv, dtype=float)
                    fig2.add_trace(go.Scatter3d(
                        x=x_arr, y=y_arr, z=z_arr,
                        mode='markers',
                        marker=dict(size=marker_size, color=norm,
                                    colorscale='Viridis', cmin=0, cmax=1, opacity=0.9,
                                    colorbar=dict(title=f"{gene} expr", len=0.5)),
                        name=f"{gene} expression", hoverinfo="name"
                    ))
            continue

        # 2–3 genes: RGB
        genes_rgb = [g for g in selected_genes[:3] if g in mtx.var_names]
        if genes_rgb:
            gv_by_gene = {g: _gene_vec(mtx, g) for g in genes_rgb}
            if child_subs:
                for sub in child_subs:
                    sub_mask = mask_cells & (mtx.obs['_subcl'] == sub)
                    sub_np = sub_mask.to_numpy()
                    x_sub = mtx.obs.loc[sub_mask, X_COL].to_numpy()
                    y_sub = scale_y(mtx.obs.loc[sub_mask, Y_COL].to_numpy())
                    z_sub = mtx.obs.loc[sub_mask, Z_COL].to_numpy()
                    N = x_sub.shape[0]
                    if N == 0: continue
                    rgb = np.zeros((N,3), dtype=float); alpha = np.zeros(N, dtype=float)
                    for chan, gene in enumerate(genes_rgb):
                        gv = gv_by_gene[gene][sub_np]
                        if gv.size and np.nanmax(gv) > 0:
                            gv = np.log1p(gv); gmin, gmax = float(np.nanmin(gv)), float(np.nanmax(gv))
                            norm = ((gv - gmin) / (gmax - gmin + 1e-9)) ** gamma
                        else:
                            norm = np.zeros_like(gv)
                        rgb[:,chan] = norm; alpha += norm
                    alpha = np.clip(alpha / max(1,len(genes_rgb)), 0.10, 1.0)
                    rgba = [f'rgba({int(255*r)},{int(255*g)},{int(255*b)},{a:.3f})'
                            for (r,g,b), a in zip(rgb, alpha)]
                    fig2.add_trace(go.Scatter3d(
                        x=x_sub, y=y_sub, z=z_sub, mode='markers',
                        marker=dict(size=marker_size, color=rgb), # use RGB, without alpha
                        name=f"RGB in {sub}", hoverinfo="name"
                    ))
            elif not all_child_subs:
                N = x_arr.shape[0]
                rgb = np.zeros((N,3), dtype=float); alpha = np.zeros(N, dtype=float)
                mask_np = mask_cells.to_numpy()
                for chan, gene in enumerate(genes_rgb):
                    gv = gv_by_gene[gene][mask_np]
                    if gv.size and np.nanmax(gv) > 0:
                        gv = np.log1p(gv); gmin, gmax = float(np.nanmin(gv)), float(np.nanmax(gv))
                        norm = ((gv - gmin) / (gmax - gmin + 1e-9)) ** gamma
                    else:
                        norm = np.zeros_like(gv)
                    rgb[:,chan] = norm; alpha += norm
                alpha = np.clip(alpha / max(1,len(genes_rgb)), 0.10, 1.0)
                rgba = [f'rgba({int(255*r)},{int(255*g)},{int(255*b)},{a:.3f})'
                        for (r,g,b), a in zip(rgb, alpha)]
                fig2.add_trace(go.Scatter3d(
                    x=x_arr, y=y_arr, z=z_arr, mode='markers',
                    marker=dict(size=marker_size, color=rgba),
                    name="RGB Expression", hoverinfo="name"
                ))
                
    fig2.update_layout(
        scene=dict(
            xaxis=dict(title='', showbackground=False, showticklabels=False),
            yaxis=dict(title='', showbackground=False, showticklabels=False),
            zaxis=dict(title='', showbackground=False, showticklabels=False),
            aspectmode='data'   # keep true units
        ),
        margin=dict(l=0,r=0,b=0,t=30),
        showlegend=False,
        uirevision='constant'
    )
    #if camera_state:
    #    fig2.update_layout(scene_camera=camera_state)


    # RGB legend
    # if len(selected_genes) > 1:
    #     color_names = ['Red','Green','Blue']
    #     items = []
    #     for i, gene in enumerate(selected_genes[:3]):
    #         items.append(html.Span([
    #             html.Span(style={'display':'inline-block','width':'16px','height':'16px',
    #                              'backgroundColor': color_names[i].lower(),'marginRight':'8px'}),
    #             f"{color_names[i]}: {gene}"
    #         ], style={'marginRight':'16px'}))
    #     rgb_legend = html.Div([html.Strong("RGB mapping: "), *items])
    # else:
    #     rgb_legend = ""
    if len(selected_genes) > 1:
        genes_rgb = selected_genes[:3]
        colors = ['red', 'green', 'blue']
        
        positions = [
            {'left': '0px',  'top': '10px'},   # red
            {'left': '22px', 'top': '10px'},   # green
            {'left': '11px', 'top': '-4px'},   # blue (above, centered)
        ]

        circles = []
        for i, gene in enumerate(genes_rgb):
            circles.append(
                html.Div(
                    gene,
                    style={
                        'position': 'absolute',
                        'left': positions[i]['left'],
                        'top': positions[i]['top'],
                        'width': '36px',
                        'height': '36px',
                        'borderRadius': '50%',
                        'backgroundColor': colors[i],
                        'opacity': '0.55',
                        'display': 'flex',
                        'alignItems': 'center',
                        'justifyContent': 'center',
                        'fontSize': '8px',
                        'color': 'white',
                        'fontWeight': 'bold',
                        'textAlign': 'center',
                        'lineHeight': '1.1',
                        'padding': '2px',
                        'boxSizing': 'border-box',
                    }
                )
            )

        # Labels next to the venn
        labels = []
        for i, gene in enumerate(genes_rgb):
            labels.append(
                html.Div([
                    html.Span(style={
                        'display': 'inline-block',
                        'width': '10px', 'height': '10px',
                        'borderRadius': '50%',
                        'backgroundColor': colors[i],
                        'marginRight': '5px',
                        'verticalAlign': 'middle'
                    }),
                    html.Span(gene, style={'verticalAlign': 'middle'})
                ], style={'marginBottom': '3px', 'fontSize': '12px'})
            )

        rgb_legend = html.Div([
            html.Div(
                circles,
                style={
                    'position': 'relative',
                    'width': '70px',
                    'height': '60px',
                    'display': 'inline-block',
                    'verticalAlign': 'middle',
                    'marginRight': '16px'
                }
            ),
            html.Div(labels, style={
                'display': 'inline-block',
                'verticalAlign': 'middle'
            })
        ], style={
            'display': 'flex',
            'alignItems': 'center',
            'justifyContent': 'center',
            'paddingTop': '8px'
        })
    else:
        rgb_legend = ""
    if gene_only_change:
        return dash.no_update, fig2, rgb_legend
    return fig1, fig2, rgb_legend



@app.callback(
    Output('section_panels', 'children'),
    Output('section_panels', 'style'),
    Input('tabs', 'value'),
    Input('cluster_selector', 'data'),
    Input('section_selector_t2', 'value'),
    Input('size_slider_t2', 'value'),
    Input('gene_selector', 'value'),
    Input('gamma_slider', 'value'),
    Input({'type': 'subcluster-check', 'index': ALL}, 'value'),
    Input('cluster_colors_store', 'data'),
    Input('subcluster_colors_store', 'data'),
    State({'type': 'subcluster-check', 'index': ALL}, 'id'),
    State({'type': 'cluster-check', 'index': ALL}, 'id'),
)
def update_section_panels(active_tab, selected_clusters, selected_sections,
                          marker_size, selected_genes, gamma,
                          sub_values, parent_colors, sub_colors,
                          sub_ids, cluster_ids):
    if active_tab != 'tab-2':
        raise PreventUpdate

    selected_clusters = selected_clusters or []
    selected_sections = [str(s) for s in (selected_sections or [])]
    selected_genes = [g for g in (selected_genes or []) if g]

    selected_subs = {
        sid['index']
        for val, sid in zip(sub_values or [], sub_ids or [])
        if (val and len(val) > 0)
    }

    live_cluster_colors = dict(parent_colors or cluster_colors)
    live_subcluster_colors = dict(sub_colors or subcluster_colors)

    gv_by_gene = {g: _gene_vec(mtx, g) for g in selected_genes if g in mtx.var_names}
    gv_norm_by_gene = {}
    for g, gv in gv_by_gene.items():
        gv_log = np.log1p(gv)
        gmin, gmax = float(np.nanmin(gv_log)), float(np.nanmax(gv_log))
        gv_norm_by_gene[g] = ((gv_log - gmin) / (gmax - gmin + 1e-9)) ** gamma

    section_cluster_idx = mtx.obs.groupby(["_section", "_cl"]).indices
    section_subcl_idx = mtx.obs.groupby(["_section", "_subcl"]).indices

    x_all = mtx.obs[X_COL].to_numpy()
    y_all = mtx.obs[Z_COL].to_numpy()

    print("section keys sample:", list(section_cluster_idx.keys())[:3])
    print("selected_sections:", selected_sections)
    print("selected_clusters:", selected_clusters)

    x_range_vals, y_range_vals = [], []
    for section in selected_sections:
        for cluster in selected_clusters:
            idx = section_cluster_idx.get((section, cluster))
            if idx is not None and len(idx):
                x_range_vals.append(x_all[idx])
                y_range_vals.append(y_all[idx])

    if x_range_vals:
        x_concat = np.concatenate(x_range_vals)
        y_concat = np.concatenate(y_range_vals)
        x_range = [float(x_concat.min()), float(x_concat.max())]
        y_range = [float(y_concat.min()), float(y_concat.max())]
    else:
        x_range, y_range = None, None

    figs = []
    for section in selected_sections:
        fig = go.Figure()

        for cluster in selected_clusters:
            idx = section_cluster_idx.get((section, cluster))
            if idx is None or len(idx) == 0:
                continue

            all_child_subs = parent_to_sub.get(cluster, [])
            child_subs = [s for s in all_child_subs if s in selected_subs]

            x_arr = x_all[idx]
            y_arr = y_all[idx]

            if len(selected_genes) == 0:
                if child_subs:
                    for sub in child_subs:
                        sub_idx = section_subcl_idx.get((section, sub))
                        if sub_idx is None or len(sub_idx) == 0:
                            continue
                        fig.add_trace(go.Scatter(
                            x=x_all[sub_idx], y=y_all[sub_idx],
                            mode='markers',
                            marker=dict(
                                size=marker_size,
                                color=live_subcluster_colors.get(sub, live_cluster_colors.get(cluster, '#808080')),
                                opacity=0.85
                            ),
                            name=f'{sub} cells', hoverinfo="name"
                        ))
                else:
                    fig.add_trace(go.Scatter(
                        x=x_arr, y=y_arr,
                        mode='markers',
                        marker=dict(
                            size=marker_size,
                            color=live_cluster_colors.get(cluster, '#808080'),
                            opacity=0.85
                        ),
                        name=f'{cluster} cells', hoverinfo="name"
                    ))
                continue

            if len(selected_genes) == 1:
                gene = selected_genes[0]
                if gene in gv_norm_by_gene:
                    if child_subs:
                        for sub in child_subs:
                            sub_idx = section_subcl_idx.get((section, sub))
                            if sub_idx is None or len(sub_idx) == 0:
                                continue
                            norm = gv_norm_by_gene[gene][sub_idx]
                            fig.add_trace(go.Scatter(
                                x=x_all[sub_idx], y=y_all[sub_idx],
                                mode='markers',
                                marker=dict(size=marker_size, color=norm,
                                            colorscale='Viridis', cmin=0, cmax=1, opacity=0.9),
                                name=f"{gene} in {sub}", hoverinfo="name"
                            ))
                    else:
                        norm = gv_norm_by_gene[gene][idx]
                        fig.add_trace(go.Scatter(
                            x=x_arr, y=y_arr,
                            mode='markers',
                            marker=dict(size=marker_size, color=norm,
                                        colorscale='Viridis', cmin=0, cmax=1, opacity=0.9),
                            name=f"{gene} expression", hoverinfo="name"
                        ))
                continue

            genes_rgb = [g for g in selected_genes[:3] if g in gv_norm_by_gene]
            if genes_rgb:
                if child_subs:
                    for sub in child_subs:
                        sub_idx = section_subcl_idx.get((section, sub))
                        if sub_idx is None or len(sub_idx) == 0:
                            continue
                        N = len(sub_idx)
                        rgb = np.zeros((N, 3), dtype=float)
                        for chan, gene in enumerate(genes_rgb):
                            rgb[:, chan] = gv_norm_by_gene[gene][sub_idx]
                        fig.add_trace(go.Scatter(
                            x=x_all[sub_idx], y=y_all[sub_idx], mode='markers',
                            marker=dict(size=marker_size, color=rgb_array_to_hex(rgb)),
                            name=f"RGB in {sub}", hoverinfo="name"
                        ))
                else:
                    N = len(idx)
                    rgb = np.zeros((N, 3), dtype=float)
                    for chan, gene in enumerate(genes_rgb):
                        rgb[:, chan] = gv_norm_by_gene[gene][idx]
                    fig.add_trace(go.Scatter(
                        x=x_arr, y=y_arr, mode='markers',
                        marker=dict(size=marker_size, color=rgb_array_to_hex(rgb)),
                        name="RGB Expression", hoverinfo="name"
                    ))

        fig.update_layout(
            xaxis=dict(visible=False, range=x_range),
            yaxis=dict(visible=False, range=y_range),
            margin=dict(l=0, r=0, b=20, t=20),
            showlegend=False,
            title=dict(text=f"Section {section}", x=0.5),
            uirevision=str(section)
        )

        figs.append(fig)

    n = len(selected_sections)
    cols = min(n, 5) if n else 1
    panel_height = 200 if n > 4 else 320 if n > 2 else 480

    grid_style = {
        'display': 'grid',
        'gridTemplateColumns': f'repeat({cols}, 1fr)',
        'gap': '8px',
    }

    panels = [
        dcc.Graph(figure=fig, style={'height': f'{panel_height}px'}, config={'displayModeBar': False})
        for fig in figs
    ]
    return panels, grid_style



if __name__ == '__main__':
    port = int(os.environ.get('PORT', 8050))
    app.run(debug=False, host='0.0.0.0', port=port)
