# ifs125preview
Standalone Python Qt program to provide live preview of the idle mode measurements of an IFS125HR.

The program includes rudimentary FFT functionality and provides a smoothed (= low pass - filtered) interferogram 
to detect potential detector nonlinearity during alignment.

To adjust to your instrument, feel free to adjust the parameters in config.yaml and copy-paste the site specific section. The measurement command are best obtained by copying them from Opus: 

While running Check signal in the preferred preview mode, check menu Measure -> Optis Setup & Service -> Optic Communications and copy-paste the relevant command series starting with cmd.htm....
