#!/bin/bash

git checkout testing
git pull
git checkout dev
git pull
git checkout master
git pull
git push
