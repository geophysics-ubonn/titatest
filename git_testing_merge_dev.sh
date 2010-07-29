#!/bin/bash

git push
git checkout dev
git merge testing
git push
git checkout testing
