# MolViewSpec Section Template

## Markdown Cell

Explain:

- what each 3D view shows
- which structure source is used
- how colors should be interpreted
- what limitations or skipped views mean

## Code Cell Shape

```python
try:
    import base64
    import json
    from IPython.display import HTML, IFrame, display
    import molviewspec as mvs
    MVS_AVAILABLE = True
except ImportError as exc:
    MVS_AVAILABLE = False
    print(f"MolViewSpec unavailable; skipping 3D views: {exc}")


def show_mol_view(state, label, width=950, height=600):
    try:
        html = state.molstar_html()
        encoded = base64.b64encode(html.encode()).decode()
        display(HTML(f'<div style="font-weight:bold;margin:8px 0 4px">{label}</div>'))
        display(IFrame(src=f"data:text/html;base64,{encoded}", width=width, height=height))
    except AttributeError:
        mvsj = state.as_dict() if hasattr(state, "as_dict") else state.to_dict()
        encoded = base64.b64encode(json.dumps(mvsj).encode()).decode()
        url = f"https://molstar.org/viewer/?mvs-data={encoded}&mvs-data-format=mvsj"
        display(HTML(f'<b>{label}</b><br><a href="{url}" target="_blank">Open in Mol* viewer</a>'))


generated_views = []
skipped_views = []

if MVS_AVAILABLE:
    # Verify API, build states, render views, append generated/skipped labels.
    pass

print({"generated_views": generated_views, "skipped_views": skipped_views})
```

