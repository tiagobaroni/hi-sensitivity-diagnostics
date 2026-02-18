# options.py
# Changelog:
# - Reduced to minimal configuration layer for hi-sensitivity-diagnostics.
# - Standardized header under GPLv3.
# - Restricted to parameters required for physical noise modeling.

# =========================================================================
# PROJECT: hi-sensitivity-diagnostics
#
# DESCRIPTION:
# Configuration module defining processing parameters used in
# spectral sensitivity diagnostics and integration-time planning.
#
# This file intentionally contains only the minimal set of parameters
# required to compute theoretical noise limits and enforce physically
# consistent decision rules.
#
# AUTHOR: Tiago Henrique França Baroni
# EMAIL: t.baroni (at) gmail (dot) com
# LINK: https://github.com/tiagobaroni/hi-sensitivity-diagnostics
# LICENSE: GNU General Public License v3.0 (GPL-3.0)
# =========================================================================
#
# Copyright (C) 2026 Tiago Henrique França Baroni
#
# This file is part of hi-sensitivity-diagnostics.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# =========================================================================
#
# FUNCTIONAL ROLE:
#
# Provides:
#   • ProcessOptions dataclass
#   • Physical parameters for Radiometer Equation modeling
#   • Controlled override mechanism with validation
#
# These parameters influence:
#   - Theoretical noise estimation (σ_th)
#   - Thermal vs instrumental regime classification
#   - Integration-time scaling laws
#
# No spectral processing logic is implemented here.
# =========================================================================


from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any


@dataclass(frozen=True, slots=True)
class ProcessOptions:
    """
    Minimal processing/physics configuration required by the sensitivity diagnostics tool.

    Notes:
    - This module intentionally includes only the parameters consumed by hi-sensitivity-diagnostics.py.
    - All values are defaults and may be overridden via with_overrides().
    """
    t_sys_est: float = 100.0
    integ_time_est: float = 60.0
    n_pol: int = 1
    enbw_factor: float = 1.0

    def validate(self) -> None:
        """Validate physical parameters."""
        if not (self.t_sys_est > 0):
            raise ValueError("t_sys_est must be > 0")
        if not (self.integ_time_est > 0):
            raise ValueError("integ_time_est must be > 0")
        if self.n_pol < 1:
            raise ValueError("n_pol must be >= 1")
        if not (self.enbw_factor > 0):
            raise ValueError("enbw_factor must be > 0")


def with_overrides(base: ProcessOptions, **kwargs: Any) -> ProcessOptions:
    """
    Create a new ProcessOptions instance overriding selected fields.
    """
    new_proc = replace(base, **kwargs)
    new_proc.validate()
    return new_proc


DEFAULT_PROCESS_OPTIONS = ProcessOptions()
