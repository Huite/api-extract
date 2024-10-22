from typing import Any
from pathlib import Path
import glob
import inspect
import importlib
import uuid


def collect_examples(globpath: str, output_file: str) -> str:
    paths = glob.glob(globpath, recursive=True)
    content = "\n\n".join([Path(p).read_text(encoding="utf-8") for p in paths])
    Path(output_file).parent.mkdir(exist_ok=True, parents=True)
    with open(output_file, encoding="utf-8", mode="w") as f:
        f.write(content)
    return


def extract_from_rst(path: str) -> list[str]:
    """Ad hoc RST parsing."""
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()

    currentmodule = ""
    references = []
    do_parse = False
    for line in lines:
        # Start parsing if we come across a Sphinx directive.
        if line.startswith(".."):
            if line.startswith(".. autosummary::"):
                do_parse = True
            # Stop parsing if it's a .. note, etc.
            else:
                do_parse = False
            if "currentmodule::" in line:
                currentmodule = line.split("::")[-1].strip() + "."
            continue

        if do_parse and line.startswith("   "):
            stripped = line.strip()
            if stripped.startswith(":"):
                continue
            if stripped:  # skip empty strings
                references.append(f"{currentmodule}{stripped}")

    return references


def get_object(reference: str) -> Any:
    module_name, *remainder = reference.split(".")
    parent = importlib.import_module(module_name)
    obj = None
    while remainder:
        if parent is None:
            return None
        first = remainder.pop(0)
        try:
            obj = parent.__dict__[first]
        except KeyError:
            try:
                obj = getattr(parent, first)
            except Exception:
                pass
        parent = obj
    return obj


def skip(name: str, doc: str, obj: Any) -> bool:
    return (
        obj is None
        or name.startswith("_")
        or doc is None
        or isinstance(obj, (tuple, list, dict))
    )


def extract_class_attrs(cls: type, reference: str) -> dict[str, tuple[Any, str]]:
    collection = {}
    for name, obj in inspect.getmembers(cls):
        doc = inspect.getdoc(obj)
        if skip(name, doc, obj):
            continue
        collection[f"{reference}.{name}"] = (obj, doc)
    return collection


def collect_docstrings(
    references: list[str], class_descend: bool = False
) -> dict[str, tuple[Any, str]]:
    collection = {}
    for reference in references:
        obj = get_object(reference)
        doc = inspect.getdoc(obj)
        if skip("", doc, obj):
            continue
        collection[reference] = (obj, doc)
        if class_descend and inspect.isclass(obj):
            class_collection = extract_class_attrs(obj, reference)
            if class_collection:
                memberlist = "\n".join(
                    [f"   * {member}" for member in class_collection]
                )
                # Add a dummy UUID so that this list has a unique identifier
                # and is never skipped during rendering.
                collection[f"{reference} Class Members"] = (uuid.uuid1(), memberlist)
                collection.update(class_collection)
    return collection


def render_docstring_collection(collection: dict[str, tuple[Any, str]]) -> str:
    content = []
    seen = set()
    for name, (item, doc) in collection.items():
        # Skip inherited methods that we've already seen.
        if item in seen:
            continue
        content.extend([name, "\n", len(name) * "=", "\n", doc, "\n\n"])
        seen.add(item)
    return "".join(content)


def create_docstring_summary(
    input_file: str,
    output_file: str,
    class_descend: bool = False,
) -> None:
    references = extract_from_rst(input_file)
    collection = collect_docstrings(references, class_descend)
    content = render_docstring_collection(collection)
    Path(output_file).parent.mkdir(exist_ok=True, parents=True)
    with open(output_file, encoding="utf-8", mode="w") as f:
        f.write(content)
    return


if __name__ == "__main__":
    prefix = "summary"
    paths = [
        r"c:\src\imod-python\docs\api\visualize.rst",
        r"c:\src\imod-python\docs\api\prepare.rst",
        r"c:\src\imod-python\docs\api\util.rst",
        r"c:\src\imod-python\docs\api\metamod.rst",
        r"c:\src\imod-python\docs\api\wq.rst",
        r"c:\src\imod-python\docs\api\msw.rst",
        r"c:\src\imod-python\docs\api\select.rst",
        r"c:\src\imod-python\docs\api\evaluate.rst",
        r"c:\src\imod-python\docs\api\flow.rst",
    ]
    for path in paths:
        create_docstring_summary(
            input_file=path,
            output_file=f"{prefix}/imod-{Path(path).stem}.rst",
            class_descend=True,
        )
    
    path = r"c:\src\xugrid\docs\api.rst"
    create_docstring_summary(
        input_file=path,
        output_file=f"{prefix}/xugrid.rst",
        class_descend=True,
    )
    path = r"c:\src\pandamesh\docs\api\index.rst"
    create_docstring_summary(
        input_file=path,
        output_file=f"{prefix}/pandamesh.rst",
        class_descend=True,
    )

    collect_examples(
        globpath=r"c:\src\xugrid\examples\*.py",
        output_file="summary/xugrid-examples.py"
    )
    collect_examples(
        globpath=r"c:\src\pandamesh\examples\*.py",
        output_file="summary/pandamesh-examples.py"
    )
    collect_examples(
        globpath=r"c:\src\imod-python\examples\mf6\*.py",
        output_file="summary/imod-mf6-examples.py"
    )