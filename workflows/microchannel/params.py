#!/usr/bin/env python3
"""The parameters of one microchannel case, read from the unrolled
interface.json, and the quantities derived from them.

CaseParams is the single typed model the worker (py_microchannel.py), the
case builder (build_case.build) and the summary all read from. The
derived values -- period, channel count, half-channel mass flow, case
name -- are defined once here.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class CaseParams:
    project_name: str
    # geometry [mm] and rotation [rpm] -- the six build parameters
    a: float
    b: float
    R: float
    s: float
    l: float
    N_rot: float
    # wheel [mm], flow [kg/s] and temperatures [K]
    D_outer: float
    D_inner: float
    mdot_stream: float
    T_hot: float
    T_cold: float
    T_init: float
    # mesh + solver
    refinement: float
    max_co: float
    n_correctors: int
    end_time: float
    n_proc: int
    contact: str

    @classmethod
    def from_interface(cls, j: dict, default_name: str) -> "CaseParams":
        """From the unrolled interface.json (form label -> value). Raises
        ValueError on a wheel with no face (needs 0 <= D_i < D_o)."""
        p = cls(
            project_name=str(j.get("project name", default_name)),
            a=float(j["channel height a [mm]"]),
            b=float(j["half-wavelength b [mm]"]),
            R=float(j["crest fillet radius R [mm]"]),
            s=float(j["foil gauge s [mm]"]),
            l=float(j["channel length l [mm]"]),
            N_rot=float(j["wheel speed N_rot [rpm]"]),
            D_outer=float(j["wheel_outer_diameter D_o [mm]"]),
            D_inner=float(j["wheel_inner_diameter D_i [mm]"]),
            mdot_stream=float(j["stream mass flow rate mdot [kg/s]"]),
            T_hot=float(j["hot inlet temperature [K]"]),
            T_cold=float(j["cold inlet temperature [K]"]),
            T_init=float(j["initial temperature [K]"]),
            refinement=float(j["mesh refinement r"]),
            max_co=float(j["maxCo"]),
            n_correctors=int(float(j["nCorrectors"])),
            end_time=float(j["end time [s]"]),
            n_proc=max(1, min(128, int(float(j.get("number of processors",
                                                    1))))),
            contact=str(j.get("contact email", "")).strip(),
        )
        if not 0.0 <= p.D_inner < p.D_outer:
            raise ValueError(f"wheel diameters need 0 <= D_i < D_o, got "
                             f"D_i={p.D_inner:g}, D_o={p.D_outer:g} mm")
        return p

    @property
    def shape(self) -> dict:
        """The six build parameters, as the gates and the mesher take them."""
        return dict(a=self.a, b=self.b, R=self.R, s=self.s, l=self.l,
                    N_rot=self.N_rot)

    # The properties below divide by the inputs: valid only once
    # true_geometry.check_parameter_values(shape) has passed.

    @property
    def period(self) -> float:
        """Oscillation period T = 60 / N_rot [s]: one wheel revolution."""
        return 60.0 / self.N_rot

    @property
    def n_ch(self) -> int:
        """Channels on the wheel face, N_ch = pi (Do^2 - Di^2) / (4 b (a + s)),
        rounded to a whole number (at least 1): the face area
        pi (Do^2 - Di^2)/4 over b (a + s) per channel (one pitch 2b x (a + s)
        holds two channels; b is the HALF-wavelength). All lengths in mm
        (the units cancel)."""
        return max(1, round(math.pi * (self.D_outer ** 2 - self.D_inner ** 2)
                            / (4.0 * self.b * (self.a + self.s))))

    @property
    def mdot(self) -> float:
        """mdotIn of the meshed domain [kg/s]. Each stream passes half the
        wheel face, so it splits over N_ch/2 channels, and the domain is
        HALF a channel (symmetryPlane): mdotIn = mdot_stream / N_ch."""
        return self.mdot_stream / self.n_ch

    @property
    def case_name(self) -> str:
        return "case_a%g_b%g_R%g_s%g_l%g_N%g" % (
            self.a, self.b, self.R, self.s, self.l, self.N_rot)
