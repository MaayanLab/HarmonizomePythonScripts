import json
import pandas as pd
from typing import Dict, Any, List

deepdive_summaries = pd.read_pickle('deepdive_gpt4o_markdown_descriptions.pkl')

def format_deepdive(gene: str):

    summary: Dict = deepdive_summaries.get(gene, {})

    def md(props: Dict[str, Any]) -> str:
        fn = COMPONENTS.get(props.get("type"))
        return fn(props) if fn else ""

    def component_recurse(props: Dict[str, Any]) -> str:
        return "".join(md(child) for child in props.get("children", []))

    COMPONENTS = {
        "root": lambda props: component_recurse(props).strip(),

        "p": lambda props: component_recurse(props).strip() + "\n\n",

        "t": lambda props: props.get("text", ""),

        # flatten string
        "b": lambda props: component_recurse(props),
        "i": lambda props: component_recurse(props),

        # keep link text
        "a": lambda props: component_recurse(props),

        # single reference
        "f": lambda props: f"[{props.get('ref')}]",

        # grouped references
        "fg": lambda props: (
            "[" +
            ",".join(md(child) for child in props.get("children", [])) +
            "]"
        ),

        # individual ref inside fg
        "fg_f": lambda props: props.get("ref"),

        # reference range
        "fg_fs": lambda props: (
            f"{props.get('start_ref')}–{props.get('end_ref')}"
        ),

        # references group
        "rg": lambda props: (
            "\nReferences\n\n"
            #+ "-" * 10 + "\n"
            + component_recurse(props).strip()
        ),

        # reference
        "r": lambda props: (
            f"[{props.get('ref')}] "
            f"{component_recurse(props).strip()}\n"
        ),
    }

    rendered = md(summary).rstrip().split('\n\n\n', 1)[0]

    ref_to_doi = {}

    def extract_dois(node: Dict[str, Any]):
        if not isinstance(node, dict):
            return

        node_type = node.get("type")

        if node_type == "r":
            ref_number = node.get("ref")
            doi = None

            children = node.get("children", [])
            for i, child in enumerate(children):
                # look for "DOI: " text node
                if (
                    child.get("type") == "t"
                    and "DOI" in child.get("text", "")
                ):
                    if i + 1 < len(children):
                        next_child = children[i + 1]
                        if next_child.get("type") == "a":
                            doi = "".join(
                                c.get("text", "")
                                for c in next_child.get("children", [])
                                if c.get("type") == "t"
                            )
                    break

            if ref_number is not None and doi:
                ref_to_doi[int(ref_number)] = doi

        for child in node.get("children", []):
            extract_dois(child)

    extract_dois(summary)

    return rendered,ref_to_doi


def load_gmt(file, dir=''):
    gmt = {}
    with open(file, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if len(line.strip().split('\t'))==1: continue
            term, _, *geneset = line.strip().split('\t')
            if geneset:
                gmt[f"{term}{dir}"] = set(geneset)
    return gmt

def load_rummageo_gmt(file):
    gmt = {}
    with open(file, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if len(line.strip().split('\t'))==1: continue
            term, _, *geneset = line.strip().split('\t')
            if geneset:
                gmt[term] = set(geneset)
    return gmt