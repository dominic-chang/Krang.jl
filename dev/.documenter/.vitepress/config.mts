import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/Krang.jl/',// TODO: replace this in makedocs!
  title: 'Krang.jl',
  description: "A VitePress Site",
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../final_site', // This is required for MarkdownVitepress to work correctly...
  head: [['link', { rel: 'icon', href: '/public/favicon.ico' }]],
  
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
      { text: 'Getting Started', link: '/getting_started' },
      { text: 'Examples', link: '/examples/coordinate_example' },
      { text: 'API', link: '/api' }
    ],

    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Getting Started', link: '/getting_started' },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Raytracing with inclination', link: '/examples/coordinate_example' },
{ text: 'Rays parameterized with Mino time', link: '/examples/mino_time_example' }]
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
      { icon: 'medium', link: 'https://medium.com/@domchang' }
    ],
    footer: {
      message: 'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a> & <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    }
  }
})
