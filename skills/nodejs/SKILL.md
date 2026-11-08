---
name: nodejs-pro-max
description: "Node.js & TypeScript development intelligence. 12 domains, 12 stacks, 300+ entries. Actions: plan, build, create, design, implement, review, fix, improve, optimize, enhance, refactor, check Node code. Domains: pattern, package, error, perf, test, async, arch, idiom, anti, tooling, module, typescript. Stacks: web-express, web-fastify, web-hono, web-nestjs, web-nextjs, web-remix, orm-prisma, orm-drizzle, auth, database, realtime, graphql. Integrations: arch_system generator."
---

# NodeJS Pro Max - Backend Generator

## When User Requests Node.js / TypeScript Work

Follow this workflow for every Node.js-related task:

### 1. Analyze Requirements
Extract from the user's request:
- **Project type**: REST API, GraphQL, microservice, serverless, full-stack framework
- **Domain**: web, database, realtime, auth
- **Scale**: small project, enterprise, library

### Default Stack (always use unless user specifies otherwise)
| Component | Default | Package |
|-----------|---------|---------|
| **Language** | **TypeScript** | `typescript` (strict mode) |
| **Router** | Hono | `hono` |
| **Database** | PostgreSQL | `pg` |
| **ORM** | Drizzle ORM | `drizzle-orm` |
| **Validation** | Zod | `zod` |
| **Testing** | Vitest | `vitest` |
| **Deploy** | Docker | `Dockerfile` |

**Project structure**: Domain-driven layout:
```
src/index.ts
src/app.ts
src/config/      — env config & setup (zod validation)
src/routes/      — API route definitions
src/controllers/ — Request handlers
src/services/    — business logic
src/models/      — TS interfaces / DB schema
src/middlewares/ — auth, error, logging
src/utils/       — helpers
package.json
tsconfig.json
Dockerfile
```

### 2. Generate Architecture System (REQUIRED for new projects)
```bash
python3 skills/nodejs/scripts/search.py "<project description>" --arch-system -p "ProjectName"
```
This gives you: architecture, packages, error strategy, patterns, structure, testing approach.

### 3. Supplement with Domain Searches
Based on the project needs, search relevant domains:
```bash
# Patterns for the task
python3 skills/nodejs/scripts/search.py "repository pattern" --domain pattern

# Package recommendations
python3 skills/nodejs/scripts/search.py "schema validation" --domain package

# Async / Concurrency patterns
python3 skills/nodejs/scripts/search.py "promise all settled" --domain async

# Performance tips
python3 skills/nodejs/scripts/search.py "event loop blocking" --domain perf

# Typescript best practices
python3 skills/nodejs/scripts/search.py "utility types" --domain typescript
```

### 4. Stack Guidelines
If using a specific framework/stack:
```bash
python3 skills/nodejs/scripts/search.py "global error handler" --stack web-express
python3 skills/nodejs/scripts/search.py "server components" --stack web-nextjs
python3 skills/nodejs/scripts/search.py "relations eager loading" --stack orm-prisma
```

### 5. Implement
Apply all recommendations from the search results:
- Use recommended architecture and project structure
- Apply TS patterns (Discriminated unions, Utility types)
- Follow error handling strategy (Centralized error handler, operational vs programmer errors)
- Implement proper async patterns (avoid `.then()`, use `async/await`)
- Add tests using Vitest or Jest.

### 6. Always Generate Docker Files
Every project MUST include these files:

**Dockerfile** (multi-stage build):
```dockerfile
# Builder stage
FROM node:20-alpine AS builder
WORKDIR /app
COPY package*.json ./
RUN npm ci
COPY . .
RUN npm run build

# Runtime stage
FROM node:20-alpine
WORKDIR /app
COPY package*.json ./
RUN npm ci --omit=dev
COPY --from=builder /app/dist ./dist
EXPOSE 3000
CMD ["node", "dist/index.js"]
```

## Pre-Delivery Checklist

Before delivering Node.js/TS code, verify:
- [ ] `tsc --noEmit` passes
- [ ] `npm run lint` catches no warnings
- [ ] `npm run test` passes
- [ ] No `any` types used in new code (use `unknown` or specific types)
- [ ] Async functions properly awaited; no floating promises
- [ ] Unhandled promise rejections properly caught
- [ ] Centralized error handler catches routes errors
- [ ] `.env` variables validated on startup (e.g. using Zod)

## Quick Domain Reference

| Domain | When to Search |
|--------|---------------|
| `pattern` | Need a design pattern (repository, dependency injection) |
| `package` | Choosing a library/package for a use case |
| `error` | Error handling strategy, custom error classes |
| `perf` | Performance optimization, memory leaks, event loop |
| `test` | Testing approach, mocking, e2e vs unit |
| `async` | Promises, event emitters, streams, worker threads |
| `typescript` | TS features, generic types, strict mode |
| `arch` | Project structure, clean architecture |
| `idiom` | Node.js conventions, naming |
| `anti` | Avoiding common mistakes |
| `tooling` | ESLint, Prettier, tsconfig options |
| `module` | ES modules vs CommonJS |

## Resources

- **scripts/**: BM25 search engine, architecture system generator
- **data/**: 12 domains and 12 stacks with curated Node.js/TS entries
