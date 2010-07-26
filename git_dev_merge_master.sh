#!/bin/bash

git push
git checkout dev
git merge master
git push
git checkout testing
git merge master
git push
git checkout master
