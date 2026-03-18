#!/usr/bin/env python3
"""Generate graphical abstract using Google Gemini Imagen API."""
import google.generativeai as genai
import os
import sys

API_KEY = os.environ.get("GOOGLE_API_KEY", "AIzaSyB-au-8VkzJefKTp84E5mXO48H-es6tXM4")
genai.configure(api_key=API_KEY)

PROMPT = """A professional, clean scientific graphical abstract for a genomics research paper published in BMC Genomics.

LEFT SIDE: A stylized computer monitor showing DNA sequence analysis with colored base pairs (A,T,G,C) and a small phylogenetic tree icon. Above it, the text "SAFARI-IGHJ" in a modern sans-serif font. Below: silhouettes of 5 wild ruminants (saiga antelope, bison, oryx, waterbuck, deer) in a row.

CENTER: A large horizontal arrow flowing from left to right, labeled "Prediction → Validation". The arrow transitions from blue (computational) to teal (biological).

RIGHT SIDE: Three stacked panels representing the 3 validation pillars:
1. Top panel: RNA helix icon with "11.4M reads" text (Transcription)
2. Middle panel: Intron-exon splice diagram with "1.49M spliced" text (Splicing)
3. Bottom panel: DNA double helix with a break/junction symbol and "1,424 V(D)J junctions" text (Recombination)

BOTTOM BANNER: "13 species | 7 tribes | 30+ individuals | 3 validation pillars"

Style: Modern flat vector illustration, blue (#2e86c1) and teal (#1abc9c) color palette with orange (#e67e22) accents. White background. Professional scientific aesthetic. No photorealistic elements. Clean lines and minimal detail. Publication-quality."""

print("Attempting image generation with Gemini...")

try:
    # Try Imagen model first
    model = genai.ImageGenerationModel("imagen-3.0-generate-002")
    result = model.generate_images(
        prompt=PROMPT,
        number_of_images=1,
        aspect_ratio="16:9",
    )
    result.images[0].save("../figures/graphical_abstract.png")
    print("SUCCESS: graphical_abstract.png saved!")
except Exception as e1:
    print(f"Imagen failed: {e1}")
    try:
        # Fallback: try gemini-2.0-flash with image generation
        model = genai.GenerativeModel("gemini-2.0-flash-exp-image-generation")
        response = model.generate_content(
            PROMPT,
            generation_config=genai.GenerationConfig(response_mime_type="image/png")
        )
        if response.candidates and response.candidates[0].content.parts:
            for part in response.candidates[0].content.parts:
                if hasattr(part, 'inline_data') and part.inline_data:
                    with open("../figures/graphical_abstract.png", "wb") as f:
                        f.write(part.inline_data.data)
                    print("SUCCESS: graphical_abstract.png saved via Gemini Flash!")
                    sys.exit(0)
        print("No image in response. Trying text-based Gemini...")
    except Exception as e2:
        print(f"Gemini Flash failed: {e2}")

    try:
        # Last fallback: regular Gemini for SVG description
        model = genai.GenerativeModel("gemini-2.0-flash")
        response = model.generate_content(
            "Generate an SVG (XML) graphical abstract for a scientific paper about SAFARI-IGHJ, a bioinformatics pipeline that discovers immune genes in wild ruminant genomes. Use blue and teal colors. Include: computer icon, DNA helix, arrow from 'Prediction' to 'Validation', and text '13 species | 7 tribes | 11.4M reads'. Return ONLY the SVG XML code, nothing else."
        )
        svg_text = response.text
        if "<svg" in svg_text:
            start = svg_text.index("<svg")
            end = svg_text.rindex("</svg>") + 6
            svg_clean = svg_text[start:end]
            with open("../figures/graphical_abstract.svg", "w") as f:
                f.write(svg_clean)
            print("SUCCESS: graphical_abstract.svg saved (SVG format)")
        else:
            print("No SVG generated. Please use the prompt in graphical_abstract_prompt.txt manually.")
    except Exception as e3:
        print(f"All methods failed: {e3}")
        print("Please generate manually using graphical_abstract_prompt.txt")
