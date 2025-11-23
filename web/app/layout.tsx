import type { Metadata } from "next";
import { Geist, Geist_Mono } from "next/font/google";
import "./globals.css";

const geistSans = Geist({
  variable: "--font-geist-sans",
  subsets: ["latin"],
});

const geistMono = Geist_Mono({
  variable: "--font-geist-mono",
  subsets: ["latin"],
});

export const metadata: Metadata = {
  title: "MoleculeAI - Quantum-Enhanced Drug Discovery",
  description: "AI-powered molecular analysis platform combining classical machine learning and quantum computing for accelerated drug discovery. Analyze molecules, predict properties, and discover new therapeutic compounds.",
  keywords: ["drug discovery", "molecular analysis", "quantum computing", "AI", "machine learning", "chemistry", "pharmaceuticals"],
  authors: [{ name: "MoleculeAI" }],
  openGraph: {
    title: "MoleculeAI - Quantum-Enhanced Drug Discovery",
    description: "AI-powered molecular analysis platform combining classical machine learning and quantum computing for accelerated drug discovery.",
    type: "website",
    siteName: "MoleculeAI",
    images: [
      {
        url: "/og-image.png",
        width: 1200,
        height: 630,
        alt: "MoleculeAI - Quantum-Enhanced Drug Discovery",
      },
    ],
  },
  twitter: {
    card: "summary_large_image",
    title: "MoleculeAI - Quantum-Enhanced Drug Discovery",
    description: "AI-powered molecular analysis platform combining classical machine learning and quantum computing for accelerated drug discovery.",
    images: ["/og-image.png"],
  },
  icons: {
    icon: [
      { url: "/favicon.ico", sizes: "any" },
      { url: "/favicon-32.png", sizes: "32x32", type: "image/png" },
      { url: "/favicon-16.png", sizes: "16x16", type: "image/png" },
    ],
    apple: [
      { url: "/apple-touch-icon.png", sizes: "180x180", type: "image/png" },
    ],
  },
  manifest: "/manifest.json",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body
        className={`${geistSans.variable} ${geistMono.variable} antialiased`}
      >
        {children}
      </body>
    </html>
  );
}
