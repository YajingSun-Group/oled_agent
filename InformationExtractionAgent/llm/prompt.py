
import os
def get_text_prompt(text):
    return f"""
Please carefully analyze the provided literature content and extract structured information related to OLED materials and devices. Follow the instructions strictly:

1. All extracted data must strictly follow the original text — no speculation or inferred content.

2. If a data point is not mentioned, return its value as null.

3. If the text describes multiple emitter molecules or multiple device configurations, extract each emitter-device pair individually as entries in a list.

4. For all device performance metrics, such as EQE (external quantum efficiency), brightness, lifetime, etc., extract the maximum reported value.


5. For each emitter, extract both:
   - "emitter_name_full" — the full chemical name or standard molecule label (e.g., “4CzIPN”, “2CzPN”)
   - "emitter_name_abbreviation" — the nickname, code, or label used in the paper (e.g., “compound 1”, “TADF-A”, “EML-1”)
   - If only one form is present, leave the other as null.


6. For the emission layer, support various architectures including:
  - pure emitters,
  - host-dopant systems,
  - multi-dopant blends,
  - multi-layer emission structures. Use the structured "emission_layer_details" format to represent this information clearly.

7. Units may vary across different papers (e.g., “%” or “percent”, “cd/m²” or “cd m⁻²”).

  - Always extract the unit as written in the original text.

  - Do not normalize or convert units — preserve them exactly as stated for traceability.

  - If no unit is mentioned, return null for the unit field.
  
8. Use the unit formats shown below and do not fabricate missing information.

✅ Return the extracted result in the following JSON format:

{{
  "materials": [
    {{
      "emitter_name_full": "string or null",
      "emitter_name_abbreviation": "string or null",
      "emitter_SMILES": "Emitter molecule SMILES",
      "emission_wavelength_material": {{"value": number or null, "unit": "nm"}},
      "emission_efficiency": {{"value": number or null, "unit": "%"}},
      "emission_lifetime": {{"value": number or null, "unit": "ns"}},
      "energy_levels": {{
        "HOMO": {{"value": number or null, "unit": "eV"}},
        "LUMO": {{"value": number or null, "unit": "eV"}}
      }}
    }}
  ],
  "devices": [
    {{
      "device_structure": {{
        "anode": "string or null",
        "hole_injection_layer": "string or null",
        "hole_transport_layer": "string or null",
        "emission_layer_details": {{
          "emission_layer_type": "pure" | "host-dopant" | "multi-dopant" | "multi-layer" | null,
          "pure_emitter": "string or null",
          "host": "string or null",
          "dopants": [
            {{
              "name": "string",
              "wt_percent": number or null
            }}
          ],
          "emission_layers": [
            {{
              "layer_index": number,
              "host": "string or null",
              "dopants": [ {{ "name": "string", "wt_percent": number or null }} ],
              "pure_emitter": "string or null"
            }}
          ]
        }},
        "electron_transport_layer": "string or null",
        "electron_injection_layer": "string or null",
        "cathode": "string or null"
      }},
      "device_emission_wavelength": {{"value": number or null, "unit": "nm"}},
      "device_brightness": {{"value": number or null, "unit": "cd/m²"}},
      "turn_on_voltage": {{"value": number or null, "unit": "V"}},
      "current_efficiency": {{"value": number or null, "unit": "cd/A"}},
      "power_efficiency": {{"value": number or null, "unit": "lm/W"}},
      "maximum_EQE": {{"value": number or null, "unit": "%"}},
      "device_lifetime": {{"value": number or null, "unit": "h"}}
    }}
  ]
}}

❗ Only return the JSON result. Do not include explanations or inferred information.

Below is the text content:
{text}
"""


def get_image_molecule_identify_prompt(images_path):
  num_images = len(os.listdir(images_path))
  PROMPT = f"""
I will provide you with {num_images} images from a scientific paper. The images are numbered from 1 to {num_images}.
Your task is to identify if each image contains molecules. Please respond with a JSON object where the key is "images_with_molecules" and the value is a list of image indices that contain molecules.
If no images contain molecules, return an empty list for "images_with_molecules".

For example:
{{"images_with_molecules": [1, 3]}}  
If no images contain molecules, return:
{{"images_with_molecules": []}}

Make sure to analyze each image thoroughly and return the JSON with your response for images that contain molecules.
"""
  return PROMPT

def get_image_label_prompt(images_path):
  num_images = len(os.listdir(images_path))
  return f"""
I will give you {num_images} images. Image 1 contains both S-BN and 2S-BN molecules. Images 2-{num_images} are individual molecules. Please identify which images (from 2-{num_orig_files}) contain S-BN and 2S-BN molecules.

Return your answer in JSON format:
{{
    "S-BN": <image_index>,
    "2S-BN": <image_index>
}}

Only return the JSON object, no additional text.
"""