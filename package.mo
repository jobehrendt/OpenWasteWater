within ;
package OpenWasteWater "Modelica WasteWater Library"
extends Modelica.Icons.Package;


annotation (
  Window(
    x=0.45,
    y=0.01,
    width=0.44,
    height=0.65,
    library=1,
    autolayout=1),
  Documentation(info= "<html><head></head><body><p>This package can be used to model and simulate biological municipal
wastewater treatment plants.</p>


<p>PS: SimulationTimeUnit for all models built with the WasteWater library is days [d].</p>

<p>The WasteWater package is free software; it can be redistributed and/or modified under the terms of&nbsp;</p><p>&lt;a href=\"https://creativecommons.org/licenses/by/4.0/\"&gt;Creative Commons Attribution 4.0 International License.&lt;/a&gt;.</p>

<p>The WasteWater package currently consists of the following subpackages:
</p><ul>
<li>ASM1            - Activated Sludge Model No.1  (models 13 wastewater components and 8 biological processes)</li>
<li>ASM3            - Activated Sludge Model No.3  (models 13 wastewater components and 12 biological processes)</li>
<li>Icons           - Icon definitions for wastewater treatment components</li>
<li>ADM             - Anaerobic Digster Model</li>
<li>ASM1P           - Activated Sludge Model No.1 with P precipitation  (models 20 wastewater components and 13 biological and chemical processes)</li>
<li>ADMP            - Anaerobic Digster Model, dealing with P-Precipitation </li>
</ul>
<p></p>

<h5>Main Author</h5>
<blockquote>
Joachim Behrendt<br>
Hamburg University of Technology<br>
Institut of Wastewater Management and Water Protection<br>
B-2<br>
21073 Hamburg<br>
Germany<br>
email: <a href=\"mailto:j.behrendt@tuhh.de\">j.behrendt@tuhh.de</a><br>
</blockquote>
<p>Copyright © 2020 - 2022, Joachim Behrendt</p>
</body></html>", revisions = "<html>
<h5>version 1.2.0 (2022-04-26)</h5></html>"));
end OpenWasteWater;