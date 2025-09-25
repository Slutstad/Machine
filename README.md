# 3D-visualisering av svart hål

Detta projekt bygger en fysikdriven visualisering av ett Schwarzschild-svart hål. Alla tidigare experimentella
"pingpong"- eller RL-prototyper har plockats bort; repot fokuserar nu uteslutande på den astrofysiska modellen och
renderingen. Paketet innehåller:

* Moduler för att beskriva svart håls egenskaper, gravitationell tidsdilatation och tidvattenkrafter.
* Post-newtonska beräkningar för hur fotoner kröks (gravitationslinsning) och uppskattningar av Shapiros tidsfördröjning.
* Numeriska integratorer som simulerar planetbanor med relativistisk periheliumförskjutning runt svart hålet.
* Ett visualiseringslager baserat på Matplotlib som renderar händelsehorisonten, en ackretionsskiva, ljusstrålar och planeter i 3D.

## Kom igång

Projektet kräver Python 3.10+ samt följande paket:

```bash
pip install numpy matplotlib
```

Kör en exempelsimulering och rendera scenen till en bildfil:

```bash
python run_simulation.py --mass 10 --output render.png
```

Parametern `--mass` anger svart hålets massa i solmassor. Standardvärdet motsvarar ett stjärnmassat svart hål på 10 solmassor. Kommandot genererar även diagnostikutskrifter för bland annat:

* Numerisk kontra svagfält-approximerad ljuskrökning.
* Orbitalperioder med relativistisk korrektion.
* Periheliumprecession för planeter i olika banor.

För att visa den interaktiva Matplotlib-figuren istället för att spara den till disk kan flaggan `--show` användas.

## Användning i Google Colab

Så här importerar du projektet till en Colab-notebook:

1. Öppna en ny notebook på [Google Colab](https://colab.research.google.com/).
2. (Valfritt) Montera ditt Google Drive för att kunna spara resultat:
   ```python
   from google.colab import drive
   drive.mount('/content/drive')
   ```
3. Klona repot i Colabs arbetskatalog:
   ```bash
   %cd /content
   !git clone https://github.com/<ditt-användarnamn>/Machine.git
   %cd Machine
   ```
4. Installera projektets beroenden:
   ```bash
   !pip install numpy matplotlib
   ```
5. Kör simuleringen eller modifiera koden efter behov, till exempel:
   ```bash
   !python run_simulation.py --mass 10 --output render.png
   ```

Kom ihåg att kopiera filer du vill behålla tillbaka till Google Drive (t.ex. `!cp render.png /content/drive/MyDrive/`). Colab-sessioner är temporära och rensas när de stängs.

## Kodstruktur

```
blackhole_vis/
├── __init__.py
├── physics/
│   ├── black_hole.py         # Fysikaliska egenskaper hos Schwarzschild-svarta hål
│   ├── constants.py          # Centrala naturkonstanter
│   ├── light.py              # Gravitationslinsning, Shapiros tidsfördröjning, fotonsfär
│   └── orbits.py             # Relativistiska banor och orbitalperioder
└── visualization/
    └── scene.py              # 3D-rendering av horisont, ackretionsskiva, ljus och planeter
run_simulation.py             # Kommandoradsverktyg för att skapa och rendera scenen
```

## Vidare utveckling

* Implementera Kerr-geometri för roterande svarta hål och tillhörande frame dragging.
* Lägg till spektral färgläggning av ackretionsskivan baserat på black-body-temperaturer.
* Koppla samman projektet med interaktiva gränssnitt (t.ex. `pythreejs` eller `vispy`) för realtidsrendering.
