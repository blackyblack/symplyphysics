from typing import Collection, Optional
from pathlib import Path
import re
import json

ID_PATTERN = re.compile(r"#\s*UNIQUE_LAW_ID:\s*(\d+)\s*")


def _check_name_uniqueness(names: Collection[str]) -> None:
    name_to_count: dict[str, int] = {}
    for name in names:
        name_to_count[name] = name_to_count.get(name, 0) + 1

    repeating_names = [name for name, count in name_to_count.items() if count > 1]
    if not repeating_names:
        return

    message = "Repeating names:\n  '{}'".format("',\n  '".join(repeating_names))
    raise ValueError(message)


def _load_saved_data(saved_path: Path) -> dict[str, int]:
    if not saved_path.exists():
        return {}

    with open(saved_path, "r", encoding="utf-8") as fp:
        id_to_name = json.load(fp)

    assert isinstance(id_to_name, dict)
    _check_name_uniqueness(id_to_name.values())
    name_to_id = {name: int(id_) for id_, name in id_to_name.items()}

    return name_to_id


def _read_file_and_match(
    path: Path,
    name: str,
    pattern: re.Pattern,
) -> tuple[list[str], int | None]:
    lines: list[str] = []
    match_id = None

    with open(path, "r", encoding="utf-8") as fp:
        for line in fp:
            m = pattern.match(line)

            if not m:
                lines.append(line)
                continue

            if match_id is not None:
                raise RuntimeError(f"Multiple ids in file '{name}'")

            match_id = int(m.group(1))

    return lines, match_id


def _add_empty_lines(lines: list[str], *, nmax: int) -> None:
    count_empty = 0
    for line in reversed(lines):
        if line != "\n":
            break
        count_empty += 1

    for _ in range(nmax - count_empty):
        lines.append("\n")


def _process_ids(
    name: str,
    name_to_id: dict[str, int],
    match_id: Optional[int],
    saved_id: Optional[int],
) -> int:
    if match_id is None:
        if saved_id is None:
            # The file is new and hasn't been saved before
            new_id = max(name_to_id.values(), default=0) + 1
            name_to_id[name] = new_id
            return new_id

        # The comment has been deleted from the file, but it's in the database
        return saved_id

    if saved_id is None:
        # The file has the comment, but its name isn't in the database, perhaps it
        # was renamed; therefore we remove the old file name and replace it with the
        # new one

        saved_names = [name for name, id_ in name_to_id.items() if id_ == match_id]
        for saved_name in saved_names:
            del name_to_id[saved_name]

        name_to_id[name] = match_id
    elif saved_id != match_id:
        # Both the file and the database contain the id, but they don't match

        message = f"File '{name}' should have index {saved_id}, got {match_id} instead"
        raise RuntimeError(message)

    return match_id


def _main(
    base_dir: Optional[Path] = None,
    saved_path: Optional[Path] = None,
    exclude_dir_names: Collection[str] = ("core", "docs"),
    pattern: re.Pattern = ID_PATTERN,
) -> None:
    base_dir = base_dir or Path(__file__).parent.parent
    saved_path = saved_path or base_dir.parent.joinpath("unique_law_id.json")
    exclude_dir_paths = [base_dir.joinpath(dir_) for dir_ in exclude_dir_names]

    name_to_id = _load_saved_data(saved_path)

    for path in base_dir.rglob("*.py"):
        if any(path.is_relative_to(dir_path) for dir_path in exclude_dir_paths):
            continue

        name = str(path.relative_to(base_dir))

        lines, match_id = _read_file_and_match(path, name, pattern)
        _add_empty_lines(lines, nmax=2)

        saved_id = name_to_id.get(name)
        write_id = _process_ids(name, name_to_id, match_id, saved_id)

        lines.append(f"# UNIQUE_LAW_ID: {write_id}\n")
        with open(path, "w", encoding="utf-8") as fp:
            fp.writelines(lines)

    id_to_name = {id_: name for name, id_ in name_to_id.items()}
    with open(saved_path, "w", encoding="utf-8") as fp:
        json.dump(id_to_name, fp, indent=4)


if __name__ == "__main__":
    _main()
