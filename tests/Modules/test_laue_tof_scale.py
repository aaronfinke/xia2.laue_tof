from __future__ import annotations

from xia2.Modules.Laue_TOF.laue_tof_scale import lawless_log_sections

# The two parts of a lawless log that go into xia2.laue_tof.log, cut down
LAWLESS_LOG = """\
        lambda   w(lambda)
        2.0000     0.40729
      w(lambda):
        1.2062 |                              :                *************
        0.4073 |**                            :
               +------------------------------------------------------------
               2                                                          4  lambda (A)
               (':' marks reference wavelength 3.0000 A)

Wavelength normalization curve written to LAMBDANORM.plt

<!--SUMMARY_BEGIN--> $TEXT:Result: $$ $$
Summary data for        Project: xia2.laue_tof Crystal: XTAL Dataset: FROMDIALS
Mn(I) half-set correlation CC(1/2)         0.832     0.885     0.070
$$ <!--SUMMARY_END-->
==== Writing merged data
"""


def test_lawless_log_sections(tmp_path):
    logfile = tmp_path / "lawless.log"
    logfile.write_text(LAWLESS_LOG)
    spectrum, summary = lawless_log_sections(logfile)
    assert spectrum[0].strip() == "w(lambda):"
    assert "marks reference wavelength" in spectrum[-1]
    assert len(spectrum) == 6
    assert summary[0].startswith("Summary data for")
    assert summary[-1].startswith("Mn(I) half-set correlation")


def test_lawless_log_sections_absent(tmp_path):
    logfile = tmp_path / "lawless.log"
    logfile.write_text("<!--SUMMARY_BEGIN--> unterminated\n")
    assert lawless_log_sections(logfile) == ([], [])
