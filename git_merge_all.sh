#!/bin/bash

git push
git checkout dev
git merge testing
git push
git checkout master
git merge dev
git push
git checkout testing
