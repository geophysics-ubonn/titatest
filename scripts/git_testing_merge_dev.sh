#!/bin/bash

git checkout testing
git merge dev
git checkout omp
git merge dev
git checkout dev
