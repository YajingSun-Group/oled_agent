# Tools

This folder contains the tools used in the OLED Agent.

## OmniParser

OmniParser is an open-source tool from Microsoft that parses user interface screenshots into structured and easily understandable elements. This significantly enhances the ability of visual language models to execute operations accurately in specific interface regions.

Key Features:
- Converts UI screenshots into structured data
- Enhances operational accuracy of visual language models  
- Supports precise interface element localization
- Open source and easy to integrate

In OLED Agent, we utilize OmniParser to parse and understand operating system interfaces, such as Windows and Ubuntu, and automatically perform tasks like paper downloading and file renaming.

For usage, please read the README.md file in the OmniParser folder.

Please run the OmniParser server before using the OLED Agent.
```bash
cd OmniParser/omnitool/omniparserserver
python -m omniparserserver
```

Please note that OmniParser is not included in the OLED Agent repository. You need to install it separately.


## DECIMER

DECIMER is a tool for splitting images into multiple molecular structures. It is a part of the OLED Agent.

For usage, please read the README.md file in the DECIMER folder.

We suggest installing DECIMER using pip.
```bash
pip install decimer
```


