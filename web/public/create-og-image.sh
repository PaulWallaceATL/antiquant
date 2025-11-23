#!/bin/bash
# Script to create og-image.png with proper colors from SVG

cd "$(dirname "$0")"

# Use rsvg-convert for proper color rendering (preserves gradients)
if command -v rsvg-convert &> /dev/null; then
    rsvg-convert --width=1200 --height=630 og-image.svg > og-image.png
    echo "✅ Created og-image.png with proper RGB colors using rsvg-convert"
elif command -v magick &> /dev/null; then
    # Fallback to ImageMagick (may lose some color info)
    magick og-image.svg -colorspace sRGB -resize 1200x630 og-image.png
    echo "✅ Created og-image.png using ImageMagick (fallback)"
else
    echo "❌ Error: Need rsvg-convert or ImageMagick to generate PNG"
    exit 1
fi

