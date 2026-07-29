# Locosopa

A github action for low-code software homepages.

## Example usage

```yaml
name: build-and-deploy

on:
  schedule:
    - cron: '0 18 * * *'
  push:
    branches:
      - main
  workflow_dispatch:

permissions:
  contents: read
  pages: write
  id-token: write

jobs:
  build-and-deploy:
    environment:
      name: github-pages
      url: ${{ steps.deployment.outputs.page_url }}
    defaults:
      run:
	shell: bash -el {0}
    runs-on: ubuntu-latest
    steps:
      - name: Checkout
        uses: actions/checkout@v6

      - name: locosopa
        uses: koesterlab/locosopa@main
        with:
          config: src/config.yaml
          path: build
```

Thereby, the config file (`src/config.yaml`) should look like this:

```yaml
# identify the software project
project:
  authors: Felix Wiegand, Johannes Koester
  copyright: 2024, Koesterlab
  name: Datavzrd

# specify primary color of lines and buttons
color: blue

docs:
  # specify path to sphinx docs to be rendered (optional, remove this section if you don't want this)
  src: docs
  # specify links to show in the navbar of the docs
  links:
    Homepage: https://datavzrd.github.io

header:
  links:
    # specify list of links for the header bar
    - text: Github
      url: https://github.com/datavzrd

hero:
  # specify a hero image (can just be a big logo)
  img:
    dark: logo_dark.svg
    light: logo_light.svg
  # specifiy a list of central statements about the software
  statements:
    - A zero/low-code, interactive, visual, server-free, browser-based reporting tool for tabular datasets
  # specify a list of central links
  links:
    - text: Tutorial
      url: https://datavzrd.github.io/docs/tutorial.html

# specify a small logo for the favicon and the top left
logo:
  dark: logo_dark.svg
    light: logo_light.svg

# specify the main development repo, for auto-retrieval of contributor stats
repo:
  contributors:
    skip:
    - github-actions[bot]
    - dependabot[bot]
  name: datavzrd/datavzrd

# specify the main features of the software
features:
  - title: Feature title
    url: optional url linking to details
    desc: |
      The feature description
    img: optional-path-to-image.png
    code: |
      optional code illustrating the feature

# specify whether to generate a fediwall (remove this section if not wanted)
fediwall:
  hashtags:
    # set hashtags to show (case insensitive)
    - datavzrd
  # specify additional accounts to follow
  accounts: []

```
