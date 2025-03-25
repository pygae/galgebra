# This is a Justfile. It is a file that contains commands that can be run with the `just` command.
# To install `just`, see https://github.com/casey/just

default:
  just --list

prep-uv:
  #!/usr/bin/env bash
  # Already installed, then exit
  which uv && exit 0
  # On Mac, install with Homebrew
  which brew && brew install uv
  # On *nix
  which uv || (curl -LsSf https://astral.sh/uv/install.sh | sh)
  # On Windows, see https://docs.astral.sh/uv/getting-started/installation/

prep-dt:
  #!/usr/bin/env bash
  which uv || (echo "Please install 'uv' in your preferred way, e.g. just prep-uv" && exit 1)
  uv venv --python=python3.11 --seed
  source .venv/bin/activate
  pip install -e .
  pip install -r test_requirements.txt

# Run tests during development
# depends on prep-dt
dt TESTS="*":
  #!/usr/bin/env bash
  source .venv/bin/activate
  pytest test/test_{{TESTS}}.py

lint:
  flake8 -v
