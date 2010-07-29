#!/bin/bash

git checkout dev
git merge testing
git checkout master
git merge dev
git checkout testing
git push
