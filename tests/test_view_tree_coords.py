import json
import re
from html.parser import HTMLParser
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Phylo.BaseTree import Clade, Tree

from gei_nj import view_tree_coords
from gei_nj.construct_nj import _set_tree_colors


@pytest.fixture
def colored_tree():
    tree = Tree(root=Clade(branch_length=0., name="Root", clades=[
        Clade(branch_length=1., clades=[Clade(branch_length=.2, name="001"),
                                      Clade(branch_length=.3, name="NA")]),
        Clade(branch_length=2., name="cell3"),
    ]), rooted=True)
    obs = pd.DataFrame({"group": ["B", "A", "B"]}, index=["001", "NA", "cell3"])
    _set_tree_colors(tree, SimpleNamespace(obs=obs), "group",
                     ["#123456", "#ff000080"], "rectangular", True)
    return tree


def read_payload(path, script_id="tree-data"):
    html = path.read_text(encoding="utf-8")
    match = re.search(r'<script type="application/json" id="' + script_id + r'">(.*?)</script>',
                      html, re.DOTALL)
    return json.loads(match.group(1))


@pytest.mark.parametrize("layout", ["rectangular", "triangular", "circular"])
def test_coordinates_match_ids_and_inherit_colors(colored_tree, tmp_path, layout):
    path = tmp_path / "coords.tsv"
    path.write_text("cell3\t30\t-3\nNA\t20\t-2\n001\t10\t-1\n")
    colored_tree.plot_metadata["layout"] = layout
    before = [(node.name, node.branch_length) for node in colored_tree.find_clades()]
    result = view_tree_coords(colored_tree, path, tmp_path / "plots" / "coords.html",
                              x_label="UMAP1", y_label="UMAP2")
    data = read_payload(result)
    assert data["coords"] == {"001": [10., -1.], "NA": [20., -2.], "cell3": [30., -3.]}
    assert data["tip_colors"] == colored_tree.tip_colors
    assert data["groups"] == [{"label": "B", "color": "#123456"},
                              {"label": "A", "color": "#ff000080"}]
    assert data["labels"] == ["UMAP1", "UMAP2"]
    assert data["layout"] == layout
    state = read_payload(result, "view-state")
    assert state["selected"] == [data["tree"]["id"]]
    assert state["branchLengths"] is True
    assert state["hideTipNames"] is False
    assert state["pointSize"] == 6.5
    assert not any(state[key] for key in ["flipX", "flipY", "swapXY"])
    assert before == [(node.name, node.branch_length) for node in colored_tree.find_clades()]


def test_custom_point_size(colored_tree, tmp_path):
    path = tmp_path / "coords.tsv"
    path.write_text("001\t1\t2\nNA\t3\t4\ncell3\t5\t6\n")
    result = view_tree_coords(colored_tree, path, tmp_path / "coords.html", point_size=9)
    assert read_payload(result, "view-state")["pointSize"] == 9


@pytest.mark.parametrize("text,message", [
    ("001\t1\t2\nNA\t3\t4\n", "1 missing, 0 extra"),
    ("001\t1\t2\nNA\t3\t4\ncell3\t5\t6\nextra\t7\t8\n", "0 missing, 1 extra"),
    ("001\t1\t2\nNA\t3\t4\n001\t5\t6\n", "unique"),
    ("001\t1\nNA\t3\ncell3\t5\n", "three headerless columns"),
])
def test_invalid_coordinate_ids_and_columns(colored_tree, tmp_path, text, message):
    path = tmp_path / "coords.tsv"
    path.write_text(text)
    with pytest.raises(ValueError, match=message):
        view_tree_coords(colored_tree, path, tmp_path / "coords.html")


@pytest.mark.parametrize("value", ["nan", "inf", "-inf"])
def test_nonfinite_coordinates(colored_tree, tmp_path, value):
    path = tmp_path / "coords.tsv"
    path.write_text(f"001\t1\t2\nNA\t{value}\t4\ncell3\t5\t6\n")
    with pytest.raises(ValueError, match="finite"):
        view_tree_coords(colored_tree, path, tmp_path / "coords.html")


def test_html_escapes_script_end_without_changing_data(colored_tree, tmp_path):
    name = "cell</script><script>alert(1)</script>"
    colored_tree.get_terminals()[0].name = name
    colored_tree.tip_colors[name] = colored_tree.tip_colors.pop("001")
    colored_tree.plot_metadata["group_colors"] = {"A</script>": "#123456"}
    path = tmp_path / "coords.tsv"
    path.write_text(f"{name}\t1\t2\nNA\t3\t4\ncell3\t5\t6\n")
    result = view_tree_coords(colored_tree, path, tmp_path / "coords.html",
                              x_label="X</script>", y_label="Y")
    data = read_payload(result)
    assert name in data["coords"] and name in data["tip_colors"]
    assert data["tree"]["children"][0]["children"][0]["tip_id"] == name
    assert data["groups"][0]["label"] == "A</script>"
    assert data["labels"] == ["X</script>", "Y"]
    assert result.read_text().count("</script>") == 3


def test_html_has_no_external_assets(colored_tree, tmp_path):
    path = tmp_path / "coords.tsv"
    path.write_text("001\t1\t2\nNA\t3\t4\ncell3\t5\t6\n")
    result = view_tree_coords(colored_tree, path, tmp_path / "coords.html")

    class Resources(HTMLParser):
        def handle_starttag(self, tag, attrs):
            values = dict(attrs)
            assert not (tag in ["script", "img", "iframe"] and "src" in values)
            assert not (tag == "link" and "href" in values)

    html = result.read_text(encoding="utf-8")
    Resources().feed(html)
    assert "__GEINJ_DATA__" not in html and "__GEINJ_STATE__" not in html
    assert read_payload(result)["labels"] == ["X", "Y"]
