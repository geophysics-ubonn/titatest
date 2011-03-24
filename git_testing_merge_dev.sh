#!/bin/bash

git checkout dev
git merge testing
git checkout omp
git merge testing
git checkout testing
