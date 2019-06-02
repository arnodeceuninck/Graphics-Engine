Features

Legende
DONE	De feature stond in de opgave, is geimplementeerd en lijkt te werken.
SKIP	De feature stond in de opgave, maar is niet geimplementeerd.
XTRA	De feature stond niet in de opgave, maar is wel geimplementeerd.

DONE	1	2D L-Systemen
DONE	0.125	L-Systemen met haakjes
SKIP	0.125	L-Systemen met stochastische replacement rules
DONE	1	3D Lijntekeningen
DONE	0.75	3D Lichamen
DONE	0.25	3D L-Systemen
DONE	1	Z-Buffering met lijnen
DONE	1.5	Z-Buffering met driehoeken
DONE	0.75	3D Fractalen
DONE	0.125	BuckyBall
DONE	0.125	Menger Spons
DONE	0.3	Ambient Licht
DONE	0.4	Diffuus Licht: lichtbron op oneindig
DONE	0.4	Diffuus Licht: puntbronnen
DONE	0.4	Speculair (glanzend) licht
DONE	1	Schaduwen
SKIP	1	Texture Mapping
DONE	0.75	3D Lijntekeningen met Bollen en Cilinders
XTRA		Texture Background
XTRA		ZBuffer Visualisation


Texture Background
Voeg volgende lijn toe aan de [General]-sectie:
background = "image.bmp"
waarbij image.bmp de afbeelding is die je als achtergrond wil hebben.
Er staat een voorbeeld .ini file in de textutes folder.
Dit is een simpele feauture die ik heb toegevoegd om toch nog iets met textures te hebben.


ZBuffer Visualisation
Dit is vooral een feauture die soms handig was om te debuggen. 
Voeg volgende lijntoe aan de [General]-sectie:
ZBufferOutput = TRUE
Dit zorgt ervoor dat er voor elke lamp i een bestand LightZbuffer[i].bmp wordt aangemaakt die de ZBuffer vanuit de postie van de lamp voorsteld. De ZBuffer van de afbeelding zelf wordt ook gevisualiseerd onder de naam ZBuffer.bmp.
In de afbeeldingen is het hoe donkerder de pixel, hoe dichterbij de waarde in de ZBuffer bij het oogpunt ligt. (Merk op: de achtergrond is zwart omdat deze als allerdichtste zijn ingesteld)
Het voorbeeld voor Texture Background bevat ook de gevisualiseerde ZBuffers.
