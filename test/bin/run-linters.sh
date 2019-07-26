#!/bin/bash

mypy vartable/*.py
mypy test/*.py
mypy ./*.py

pylint --diable=C --disable=bad-indentation vartable test

test/bin/shellcheck *.sh
test/bin/shellcheck test/bin/*.sh
