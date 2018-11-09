#!/usr/bin/env bash

echo "running twine..."
twine upload -u masakistan -p $1 wheelhouse/*
