// Script to generate logo assets from SVG
// This would typically use a library like sharp or canvas, but for now we'll create a simple version
// For production, use an image processing library or convert manually

console.log(`
Logo generation instructions:

1. Use the SVG file at: web/public/logo.svg
2. Convert to PNG using an online tool or image editor:
   - favicon.ico: 32x32 or 16x16
   - apple-touch-icon.png: 180x180
   - logo-192.png: 192x192
   - logo-512.png: 512x512
   - og-image.png: 1200x630 (for social sharing)

Or use ImageMagick:
  convert -background none -resize 180x180 web/public/logo.svg web/public/apple-touch-icon.png
  convert -background none -resize 192x192 web/public/logo.svg web/public/logo-192.png
  convert -background none -resize 512x512 web/public/logo.svg web/public/logo-512.png
  convert -background '#2563eb' -resize 1200x630 -gravity center -extent 1200x630 web/public/logo.svg web/public/og-image.png

For now, creating simple placeholder files...
`);

// Note: In production, use actual image conversion
// For now, the SVG will work for favicon and the metadata will reference it

