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
from dash_auth import BasicAuth
import boto3

if not os.path.exists("../mockdata.h5ad"):
    s3 = boto3.client('s3',
        endpoint_url=os.environ['https://t3.storageapi.dev'],
        aws_access_key_id=os.environ['tid_jAZXdNQUzMgAGKJBdyxpUfsMAwyvVfnDDPdvjQzJHHmtwciGFL'],
        aws_secret_access_key=os.environ['tsec_OGI5frIs1H3ylTnsl2Wt7rU6m+VK+OaX5fpOoRN1-Nlj1fO0j1PE4D3PG1RUhZMRsfCj4g']
    )
    s3.download_file(os.environ['ample-snackbox-ozv204ht9m'], 'mockdata.h5ad', 'mockdata.h5ad')

############## define data ###############
H5AD_FILE = "mockdata.h5ad"

CLUSTER_COL = "cluster"
SUBCLUSTER_COL = "subcluster"
SECTION_COL = "section"

X_COL = "x"
Y_COL = "y"
Z_COL = "z"

COLORS_FILE = ""#"cluster_colors.txt"
COLORS_FILE_PREFIX = ""#"e"
COLORS_FILE_CLUSTER_COLUMN = ""#"cell cluster"
COLORS_FILE_COLOR_COLUMN = ""#"color code"

MESH_FOLDER = "scaffolds"

#PANEL_SPLIT = [{"label": "Embryo", "prefix": "e"},{"label": "Placenta", "prefix": "p"},]

# disable split:
PANEL_SPLIT = None



# ---------- load data ----------
mtx = sc.read_h5ad(H5AD_FILE)

cluster_column = CLUSTER_COL
subcluster_column = SUBCLUSTER_COL
section_column = SECTION_COL

section_labels = sorted(mtx.obs[section_column].unique())
gene_list = mtx.var_names.tolist()

mtx.obs["_cl"] = mtx.obs[cluster_column].astype(str)
mtx.obs["_subcl"] = mtx.obs[subcluster_column].astype(str)

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
cluster_labels = sorted(
    mtx.obs["_cl"].unique()
)

cluster_colors = {
    c: cluster_color_map.get(c, "#808080")
    for c in cluster_labels
}

cluster_groups = {}

if PANEL_SPLIT:
    for panel in PANEL_SPLIT:
        cluster_groups[panel["label"]] = sorted(
            [
                c for c in cluster_labels
                if c.startswith(panel["prefix"])
            ]
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
        parent_to_sub[parent] = sorted(
            map(str, subs)
        )
parent_map = {sub: parent for parent, subs in parent_to_sub.items() for sub in subs}
subcluster_colors = {sub: cluster_colors.get(parent_map[sub], '#808080') for sub in parent_map}

def build_cluster_panel(label, clusters):
    print(label, len(clusters))
    id = {"type": "toggle-all", "group": label}

    return html.Div([
        html.Div([
            html.Div(label, style={'fontWeight':'bold'}),
            html.Div("Surface",
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
        'minWidth':'180px'
    })
# ---------- include scaffolds ----------
mesh_folder = "scaffolds"
mesh_data = {}

def _cluster_key_from_fname(fname: str):
    m = re.match(r"(.+)_mesh\.json$", fname)  # "e_0_mesh.json"
    if m: return m.group(1)
    m = re.match(r"mesh_(.+)\.json$", fname)  # "mesh_e_0.json"
    if m: return m.group(1)
    return None

if os.path.isdir(mesh_folder):
    for fname in os.listdir(mesh_folder):
        if not fname.endswith(".json"): continue
        key = _cluster_key_from_fname(fname)
        if key is None: continue
        try:
            with open(os.path.join(mesh_folder, fname), "r") as f:
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
        className="mesh-only"
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
        ], className="mesh-row")

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
                ], className="sub-row")
            )
        blocks.append(html.Div([prow, html.Div(srows, style={'paddingLeft':'20px'})]))
    return blocks

# ---------- Dash app ----------
app = dash.Dash(__name__)
USERS = {
    os.environ.get('DASH_USER', 'user'): os.environ.get('DASH_PASSWORD', 'password')
}
auth = BasicAuth(app, USERS)

app.layout = html.Div([

    html.H1(
        "3D Cluster Viewer",
        style={'textAlign': 'center'}
    ),

    # --------------------------------------------------
    # Settings row
    # --------------------------------------------------
    html.Div([

        html.Div([
            html.Label("Gamma (GE Color Intensity):",
                       style={'fontSize': '12px'}),
            dcc.Slider(
                id='gamma_slider',
                min=0.1,
                max=2.0,
                step=0.1,
                value=1.0,
                tooltip={"placement": "bottom",
                         "always_visible": False}
            ),

            html.Label("Marker Size:",
                       style={'fontSize': '12px'}),
            dcc.Slider(
                id='size_slider',
                min=0.1,
                max=5,
                step=0.1,
                value=0.3,
                tooltip={"placement": "bottom",
                         "always_visible": False}
            ),

            html.Label("Mesh Opacity:",
                       style={'fontSize': '12px'}),
            dcc.Slider(
                id='mesh_opacity',
                min=0.05,
                max=1.0,
                step=0.05,
                value=0.25,
                tooltip={"placement": "bottom",
                         "always_visible": False}
            ),

            html.Label("Section spacing (Z Zoom):",
                       style={'fontSize': '12px'}),
            dcc.Slider(
                id='z_zoom_slider',
                min=0.25,
                max=2.0,
                step=0.05,
                value=1.0,
                tooltip={"placement": "bottom",
                         "always_visible": False}
            )
        ], style={
            'flex': '1',
            'minWidth': '220px'
        }),

        html.Div([
            html.Label("Select Section to Display:"),

            dcc.Checklist(
                id='section_selector',
                options=[
                    {'label': str(sec), 'value': sec}
                    for sec in section_labels
                ],
                value=section_labels,
                inline=True,
                style={'margin-bottom': '10px'}
            ),

            html.Label(
                "Select up to 3 Genes (for RGB coloring):",
                style={'fontWeight': 'bold'}
            ),

            dcc.Dropdown(
                id='gene_selector',
                options=[
                    {'label': g, 'value': g}
                    for g in gene_list
                ],
                multi=True,
                value=[],
                placeholder="Pick genes to color cells",
                style={'width': '100%'}
            )
        ], style={
            'flex': '1',
            'minWidth': '280px',
            'padding': '0 20px'
        })

    ], style={
        'display': 'flex',
        'gap': '20px',
        'margin-bottom': '16px',
        'alignItems': 'flex-start'
    }),

    # --------------------------------------------------
    # Main content
    # --------------------------------------------------
    html.Div([

        # graphs
        html.Div([

            html.Div([
                dcc.Graph(
                    id='cluster_3d_plot',
                    style={
                        'height': '600px',
                        'width': '100%'
                    }
                )
            ], style={
                'flex': '1',
                'minWidth': '0',
                'marginRight': '8px'
            }),

            html.Div([
                dcc.Graph(
                    id='gene_expression_plot',
                    style={
                        'height': '600px',
                        'width': '100%'
                    }
                ),

                html.Div(
                    id='rgb_legend',
                    style={
                        'textAlign': 'center',
                        'paddingTop': '6px'
                    }
                )

            ], style={
                'flex': '1',
                'minWidth': '0',
                'marginLeft': '8px'
            })

        ], style={
            'flex': '3',
            'display': 'flex',
            'flexDirection': 'row',
            'minWidth': '0'
        }),

        # cluster controls
        html.Div(
            [
                build_cluster_panel(label, clusters)
                for label, clusters in cluster_groups.items()
            ],
            style={
                'flex': '1',
                'display': 'flex',
                'flexDirection': 'row',
                'maxHeight': '600px',
                'overflowY': 'auto',
                'paddingLeft': '4px'
            }
        )

    ], style={
        'display': 'flex',
        'gap': '4px',
        'alignItems': 'flex-start'
    }),

    dcc.Store(
        id='cluster_selector',
        data=cluster_labels
    ),

    dcc.Store(
        id='cluster_colors_store',
        data=cluster_colors
    ),

    dcc.Store(
        id='subcluster_colors_store',
        data=subcluster_colors
    ),

    dcc.Store(id='camera-store')

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
@app.callback(
    Output('cluster_selector', 'data'),
    Input({'type': 'cluster-check', 'index': ALL}, 'value'),
)
def sync_selected_clusters(values):
    selected_clusters = [val for v in (values or []) for val in (v or [])]
    return selected_clusters

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
    Output('camera-store', 'data'),
    Input('cluster_3d_plot', 'relayoutData'),
    Input('gene_expression_plot', 'relayoutData'),
    prevent_initial_call=True
)
def save_camera_state(relayout_left, relayout_right):

    trigger = ctx.triggered_id

    if trigger == 'cluster_3d_plot' and relayout_left:
        return relayout_left

    if trigger == 'gene_expression_plot' and relayout_right:
        return relayout_right

    raise PreventUpdate
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

# ----------------------------------------

def _gene_vec(adata, gene):
    """Return 1-D numpy array of expression for 'gene' across all cells."""
    Xg = adata[:, gene].X
    if hasattr(Xg, "toarray"):
        arr = Xg.toarray()
    else:
        arr = np.asarray(Xg)
    return np.ravel(arr)



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
    Input({'type': 'cluster-color-input', 'index': ALL}, 'value'),
    Input({'type': 'subcluster-color-input', 'index': ALL}, 'value'),
    State({'type': 'mesh-toggle', 'index': ALL}, 'id'),
    State({'type': 'subcluster-check', 'index': ALL}, 'id'),
    State({'type': 'cluster-color-input', 'index': ALL}, 'id'),
    State({'type': 'subcluster-color-input', 'index': ALL}, 'id'),
    State('camera-store', 'data')
)
def update_figures(
    selected_clusters, selected_sections, selected_genes,
    gamma, marker_size, mesh_opacity, z_zoom,
    mesh_values, sub_values, parent_colors, sub_colors,
    mesh_ids, sub_ids, parent_color_ids, sub_color_ids, camera_state
):
    selected_clusters = selected_clusters or []

    mesh_toggle_dict = {
        cid['index']: val
        for cid, val in zip(mesh_ids or [], mesh_values or [])
    }

    live_cluster_colors = dict(cluster_colors)
    if parent_colors and parent_color_ids:
        for col, cid in zip(parent_colors, parent_color_ids):
            if col and cid and 'index' in cid:
                live_cluster_colors[cid['index']] = col

    live_subcluster_colors = dict(subcluster_colors)
    if sub_colors and sub_color_ids:
        for col, cid in zip(sub_colors, sub_color_ids):
            if col and cid and 'index' in cid:
                live_subcluster_colors[cid['index']] = col

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

    y_factor = z_zoom or 1.0
    base_mask = mtx.obs[SECTION_COL].isin(selected_sections or section_labels)
    if selected_clusters:
        base_mask &= mtx.obs['_cl'].isin(selected_clusters)

    y0 = float(np.nanmedian(mtx.obs.loc[base_mask, Y_COL])) if base_mask.any() else 0.0
    def scale_y(yvals):
        return (yvals - y0) * y_factor + y0
    
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
                    legendgroup=parent
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
                            showlegend=True
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
                    showlegend=True
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
    y_factor = z_zoom or 1.0
    base_mask = mtx.obs[SECTION_COL].isin(selected_sections or section_labels)
    if selected_clusters:
        base_mask &= mtx.obs['_cl'].isin(selected_clusters)
    y0 = float(np.nanmedian(mtx.obs.loc[base_mask, Y_COL])) if base_mask.any() else 0.0

    def scale_y(arr):
        arr = np.asarray(arr)
        return (arr - y0) * y_factor + y0

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
                    showlegend=False
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
                        showlegend=False
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
                        name=f'{sub} cells'
                    ))
            else:
                fig2.add_trace(go.Scatter3d(
                    x=x_arr, y=y_arr, z=z_arr,
                    mode='markers',
                    marker=dict(
                        size=marker_size,
                        color=live_cluster_colors.get(cluster, '#808080'),
                        opacity=0.85
                    ),
                    name=f'{cluster} cells'
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
                            name=f"{gene} in {sub}"
                        ))
                else:
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
                        name=f"{gene} expression"
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
                        name=f"RGB in {sub}"
                    ))
            else:
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
                    name="RGB Expression"
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
    if len(selected_genes) > 1:
        color_names = ['Red','Green','Blue']
        items = []
        for i, gene in enumerate(selected_genes[:3]):
            items.append(html.Span([
                html.Span(style={'display':'inline-block','width':'16px','height':'16px',
                                 'backgroundColor': color_names[i].lower(),'marginRight':'8px'}),
                f"{color_names[i]}: {gene}"
            ], style={'marginRight':'16px'}))
        rgb_legend = html.Div([html.Strong("RGB mapping: "), *items])
    else:
        rgb_legend = ""

    return fig1, fig2, rgb_legend


if __name__ == '__main__':
    port = int(os.environ.get('PORT', 8050))
    app.run(debug=False, host='0.0.0.0', port=port)
