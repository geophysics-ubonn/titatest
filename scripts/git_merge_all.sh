#!/bin/bash

git checkout dev
git merge testing
git checkout master
git merge testing
git checkout testing
git branch
