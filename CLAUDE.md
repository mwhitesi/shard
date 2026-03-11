# CLAUDE.md

## Environment

- **Python**: >=3.10
- **Package manager**: uv
- **Setup**: `uv sync` (creates `.venv` and installs all deps including dev group)
- **Add deps**: `uv add <package>` or `uv add --dev <package>`

## Commands

- **Tests**: `uv run pytest`
- **Lint**: `uv run ruff check shard/ tests/`
- **Format**: `uv run ruff format shard/ tests/`
