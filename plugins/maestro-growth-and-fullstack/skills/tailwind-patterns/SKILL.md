---
name: tailwind-patterns
description: Tailwind CSS patterns for marketing websites and landing pages.
allowed-tools: Read, Glob, Grep
---

# Tailwind Patterns for Marketing

> Essential Tailwind CSS patterns for marketing sites.

---

## 1. Hero Section Patterns

### Centered Hero

```html
<section class="relative min-h-screen flex items-center justify-center bg-gradient-to-br from-purple-900 to-indigo-900">
  <div class="text-center max-w-4xl px-6">
    <h1 class="text-5xl md:text-7xl font-bold text-white mb-6">Your Headline</h1>
    <p class="text-xl text-gray-300 mb-8 max-w-2xl mx-auto">Supporting text here</p>
    <div class="flex flex-col sm:flex-row gap-4 justify-center">
      <a href="#" class="px-8 py-4 bg-white text-purple-900 font-semibold rounded-full hover:bg-gray-100 transition">Get Started</a>
      <a href="#" class="px-8 py-4 border border-white text-white rounded-full hover:bg-white/10 transition">Learn More</a>
    </div>
  </div>
</section>
```

---

## 2. CTA Button Patterns

### Primary Button

```html
<button class="px-8 py-4 bg-gradient-to-r from-purple-600 to-indigo-600 text-white font-semibold rounded-lg shadow-lg hover:shadow-xl hover:scale-105 transition-all duration-200">
  Get Started Free
</button>
```

### Ghost Button

```html
<button class="px-8 py-4 border-2 border-gray-800 text-gray-800 font-semibold rounded-lg hover:bg-gray-800 hover:text-white transition-all duration-200">
  Learn More
</button>
```

---

## 3. Feature Grid

```html
<div class="grid md:grid-cols-3 gap-8">
  <div class="p-8 bg-white rounded-2xl shadow-lg hover:shadow-xl transition">
    <div class="w-12 h-12 bg-purple-100 rounded-lg flex items-center justify-center mb-6">
      <!-- Icon -->
    </div>
    <h3 class="text-xl font-bold mb-3">Feature Title</h3>
    <p class="text-gray-600">Feature description goes here.</p>
  </div>
</div>
```

---

## 4. Testimonial Card

```html
<div class="p-8 bg-white rounded-2xl shadow-lg">
  <div class="flex items-center gap-1 mb-4">
    <!-- 5 stars -->
    <span class="text-yellow-400">★★★★★</span>
  </div>
  <p class="text-gray-700 mb-6 italic">"Testimonial quote here..."</p>
  <div class="flex items-center gap-4">
    <img src="avatar.jpg" class="w-12 h-12 rounded-full" alt="">
    <div>
      <p class="font-semibold">John Doe</p>
      <p class="text-sm text-gray-500">CEO, Company</p>
    </div>
  </div>
</div>
```

---

## 5. Pricing Card

```html
<div class="p-8 bg-white rounded-2xl shadow-lg border-2 border-purple-600 relative">
  <span class="absolute -top-3 left-1/2 -translate-x-1/2 px-4 py-1 bg-purple-600 text-white text-sm rounded-full">
    Popular
  </span>
  <h3 class="text-2xl font-bold mb-2">Pro Plan</h3>
  <p class="text-gray-600 mb-6">For growing teams</p>
  <div class="mb-6">
    <span class="text-5xl font-bold">$49</span>
    <span class="text-gray-500">/month</span>
  </div>
  <ul class="space-y-3 mb-8">
    <li class="flex items-center gap-2">
      <span class="text-green-500">✓</span> Feature one
    </li>
  </ul>
  <button class="w-full py-4 bg-purple-600 text-white rounded-lg font-semibold hover:bg-purple-700 transition">
    Get Started
  </button>
</div>
```

---

## 6. FAQ Accordion

```html
<div class="border-b border-gray-200">
  <button class="w-full py-6 flex justify-between items-center text-left">
    <span class="font-semibold text-lg">Question here?</span>
    <span class="text-2xl">+</span>
  </button>
  <div class="pb-6 text-gray-600">
    Answer content here.
  </div>
</div>
```

---

## 7. Social Proof Bar

```html
<div class="py-12 bg-gray-50">
  <p class="text-center text-gray-500 mb-8">Trusted by 1,000+ companies</p>
  <div class="flex flex-wrap justify-center items-center gap-12 opacity-50">
    <!-- Logo images -->
  </div>
</div>
```

---

## 8. Newsletter Form

```html
<div class="max-w-md mx-auto">
  <form class="flex gap-2">
    <input 
      type="email" 
      placeholder="Enter your email" 
      class="flex-1 px-4 py-3 border border-gray-300 rounded-lg focus:ring-2 focus:ring-purple-600 focus:border-transparent outline-none"
    >
    <button class="px-6 py-3 bg-purple-600 text-white rounded-lg font-semibold hover:bg-purple-700 transition">
      Subscribe
    </button>
  </form>
</div>
```

---

## 9. Responsive Navigation

```html
<nav class="fixed top-0 w-full bg-white/80 backdrop-blur-lg z-50 border-b border-gray-200">
  <div class="max-w-7xl mx-auto px-6 py-4 flex justify-between items-center">
    <a href="#" class="text-2xl font-bold">Logo</a>
    <div class="hidden md:flex items-center gap-8">
      <a href="#" class="text-gray-600 hover:text-gray-900 transition">Features</a>
      <a href="#" class="text-gray-600 hover:text-gray-900 transition">Pricing</a>
      <a href="#" class="px-6 py-2 bg-purple-600 text-white rounded-lg hover:bg-purple-700 transition">Get Started</a>
    </div>
  </div>
</nav>
```

---

## 10. Animation Classes

### Subtle Hover

```html
class="hover:scale-105 transition-transform duration-200"
class="hover:shadow-xl transition-shadow duration-200"
class="hover:-translate-y-1 transition-transform duration-200"
```

### Entrance Animations

```html
class="animate-fade-in"
class="animate-slide-up"
class="animate-bounce-in"
```

---

> **Remember:** Keep utility classes organized. Group layout → spacing → typography → colors → effects.
