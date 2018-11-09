#!/usr/bin/env bash

echo "running twine..."
twine upload -u masakistan -p $0 wheelhouse/*
