#!/bin/bash
# Script to create og-image.png with proper gradient background

cd "$(dirname "$0")"

# Create gradient background
magick -size 1200x630 gradient:'rgb(37,99,235)-rgb(147,51,234)' bg.png

# Composite SVG over gradient
magick bg.png og-image.svg -gravity center -composite og-image.png

# Clean up
rm bg.png

echo "✅ Created og-image.png with gradient background"

