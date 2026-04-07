# Default target
default:
  @just --list

# Install dependencies (dev group)
install:
  uv sync

# Install without dev dependencies
install-prod:
  uv sync --no-dev

# Alias for install
sync:
  uv sync

# Run tests
test:
  uv run pytest tests/ -v

# Run tests with coverage
test-cov:
  uv run pytest tests --cov=src/tdfextractor --cov-branch --cov-report=term-missing --cov-report=html --cov-report=xml

# Lint with ruff
lint:
  uv run ruff check src/ tests/

# Format code (import sort + unused imports + ruff format)
format:
  uv run ruff check --select I --fix src/ tests/
  uv run ruff check --select F401 --fix src/ tests/
  uv run ruff format src/ tests/

# Type check with ty
ty:
  uv run ty check src/

# Run lint, type check, and tests
check:
  just lint
  just ty
  just test

# Build wheel and sdist
build:
  uv build

# Upgrade Python syntax to 3.12+
upgrade:
  @echo "Upgrading Python syntax to 3.12+..."
  -@find src/tdfextractor tests -name "*.py" -type f -exec uv run pyupgrade --py312-plus {} +
  @echo "Python syntax upgraded to 3.12+"

# Publish to PyPI (requires UV_PUBLISH_TOKEN or interactive auth)
publish: build
  uv publish

# Clean build artifacts and caches
clean:
  rm -rf build/
  rm -rf dist/
  rm -rf *.egg-info
  rm -rf src/*.egg-info
  rm -rf .pytest_cache
  rm -rf .mypy_cache
  rm -rf .ruff_cache
  rm -rf htmlcov
  rm -f coverage.xml junit.xml
  find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
  find . -type f -name "*.pyc" -delete
