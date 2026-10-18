#!/usr/bin/env python

"""Tests for `hydrotools` package."""

import pytest
from hydrotools import hydrotools


def test_point_in_poly():
    polygon = [(0, 10), (10, 10), (10, 0), (0, 0)]

    point_x = 5
    point_y = 5
    assert hydrotools.point_in_poly(point_x, point_y, polygon) is True

    assert hydrotools.point_in_poly(11, 10, polygon) is False
