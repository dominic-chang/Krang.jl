import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/Krang.jl/v0.4/',// TODO: replace this in makedocs!
  title: 'Krang.jl',
  description: "Kerr Ray-tracer for Analytic Null Geodesics",
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../2', // This is required for MarkdownVitepress to work correctly...
  head: [['link', { rel: 'icon', href: '/Krang_logo.ico' }]],
  
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    // https://vitepress.dev/reference/default-theme-config
    logo: { src: '/Krang_logo.png', width: 24, height: 24 },
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Getting Started', link: '/what_is_krang' },
      { text: 'Examples', link: '/examples/coordinate-example' },
      { text: 'API', link: '/api' }
    ],
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Getting Started', collapsed: false, items: [
{ text: 'What is Krang.jl?', link: '/what_is_krang' },
{ text: 'Meshes', link: '/meshes' },
{ text: 'Raytracing', link: '/raytracing' }]
 },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Coordinates with inclination (θs)', link: '/examples/coordinate-example' },
{ text: 'Defining a Custom Material and Saving Output', link: '/examples/custom-material-example' },
{ text: 'Raytracing a Level Set geometry', link: '/examples/level-set-example' },
{ text: 'Rays & Coordinate information with Mino time (τ)', link: '/examples/mino-time-example' },
{ text: 'Neural Network Emission Model Example', link: '/examples/neural-net-example' },
{ text: 'Creating a Custom Dual Cone Model', link: '/examples/polarization-example' },
{ text: 'Raytracing a polygon mesh', link: '/examples/raytracing-mesh-example' }]
 },
{ text: 'Theory', collapsed: false, items: [
{ text: 'Raytracing', collapsed: false, items: [
{ text: 'The Geodesics of the Kerr Metric', link: '/kerr_geodesic_summary' },
{ text: 'Regularized time', link: '/time_regularization' }]
 },
{ text: 'Polarization', collapsed: false, items: [
{ text: 'Newmann-Penrose Formalism', link: '/newmann_penrose' },
{ text: 'Linear-Polarization', link: '/polarization' }]
 }]
 },
{ text: 'Krang api', link: '/api' }
]
,
    editLink: {
      pattern: 'https://github.com/LuxDL/DocumenterVitepress.jl/edit/master/docs/src/:path' // TODO: replace this in makedocs!
    },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/dchang10/Krang.jl' },
      { icon: {
        svg: '<svg xmlns="http://www.w3.org/2000/svg" x="0px" y="0px" width="48" height="48" viewBox="0 0 48 48" style="fill:#000000;"> <circle cx="14" cy="24" r="12" fill="#717171"></circle><ellipse cx="34" cy="24" fill="#717171" rx="6" ry="11"></ellipse><ellipse cx="44" cy="24" fill="#717171" rx="2" ry="10"></ellipse></svg>'
      }, link: 'https://medium.com/@domchang' }
    ],
    footer: {
      message: 'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a> & <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    }
  }
})
