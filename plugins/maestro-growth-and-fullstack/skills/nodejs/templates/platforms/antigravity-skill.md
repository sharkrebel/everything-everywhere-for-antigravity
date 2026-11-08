---
name: NodeSkill
description: AI agent that builds production-ready Node.js/TypeScript backends from natural language. Search Node.js best practices, analyze project ideas, generate complete projects — all locally.
---

# NodeSkill — Node.js/TypeScript Backend Generator

## When User Requests Node.js/TypeScript Work

Follow this workflow for every Node.js/TypeScript-related task:

### 1. Analyze Requirements
Extract from the user's request:
- **Project type**: REST API, full-stack, microservice, CLI, serverless
- **Domain**: web, database, GraphQL, real-time, auth
- **Constraints**: specific framework? serverless? edge runtime?
- **Scale**: small project, enterprise, library

### Default Stack (always use unless user specifies otherwise)
| Component | Default | Package |
|-----------|---------|---------|
| **Language** | TypeScript 5.x | `typescript` |
| **Runtime** | Node.js 20+ | `node` |
| **Web Framework** | Express / Fastify | `express`, `fastify` |
| **ORM** | Prisma | `@prisma/client` |
| **Validation** | Zod | `zod` |
| **Auth** | JWT | `jsonwebtoken` |
| **Logging** | Pino | `pino`, `pino-http` |
| **Testing** | Vitest | `vitest` |
| **Deploy** | Docker | `Dockerfile` |

**Project structure**: Modular layout:
```
src/
├── index.ts           — Entry point, server bootstrap
├── config/env.ts      — Zod-validated environment variables
├── middleware/         — Auth, error handler, logging
├── lib/               — Shared utilities (prisma client, logger)
└── modules/{entity}/
    ├── {entity}.schema.ts      — Zod schemas
    ├── {entity}.controller.ts  — Route handlers
    ├── {entity}.service.ts     — Business logic
    └── {entity}.repository.ts  — Database access (Prisma)
prisma/schema.prisma   — Database schema
```

### 2. Search Best Practices (REQUIRED)
```bash
# Search relevant domains for the task
python3 .agent/skills/nodeskill/scripts/search.py "<description>" --domain pattern
```

### 3. Supplement with Domain Searches
Based on the project needs, search relevant domains:
```bash
# Design patterns
python3 .agent/skills/nodeskill/scripts/search.py "middleware factory" --domain pattern

# Package recommendations
python3 .agent/skills/nodeskill/scripts/search.py "authentication library" --domain package

# Error handling
python3 .agent/skills/nodeskill/scripts/search.py "custom error class" --domain error

# Async patterns
python3 .agent/skills/nodeskill/scripts/search.py "promise parallel streams" --domain async

# TypeScript type patterns
python3 .agent/skills/nodeskill/scripts/search.py "discriminated union branded" --domain typescript

# Performance optimization
python3 .agent/skills/nodeskill/scripts/search.py "event loop memory" --domain perf

# Testing patterns
python3 .agent/skills/nodeskill/scripts/search.py "integration vitest mock" --domain test
```

### 4. Stack Guidelines
If using a specific framework/stack:
```bash
python3 .agent/skills/nodeskill/scripts/search.py "middleware error handling" --stack web-express
python3 .agent/skills/nodeskill/scripts/search.py "plugin schema validation" --stack web-fastify
python3 .agent/skills/nodeskill/scripts/search.py "guards interceptors pipes" --stack web-nestjs
python3 .agent/skills/nodeskill/scripts/search.py "server components actions" --stack web-nextjs
python3 .agent/skills/nodeskill/scripts/search.py "relations transactions" --stack orm-prisma
python3 .agent/skills/nodeskill/scripts/search.py "schema queries joins" --stack orm-drizzle
python3 .agent/skills/nodeskill/scripts/search.py "dataloader resolvers" --stack graphql
python3 .agent/skills/nodeskill/scripts/search.py "socket.io rooms events" --stack realtime
```

### 5. Implement
Apply all recommendations from the search results:
- Use recommended architecture and project structure
- Apply patterns (Factory, DI, Middleware Chain, etc.)
- Follow error handling strategy (Custom Error classes, Result pattern)
- Implement proper async patterns (Promise.all, Streams, AbortController)
- Use TypeScript strict mode with proper type patterns
- Add tests following testing strategy (Vitest, supertest)

### 6. Always Generate Docker Files
Every project MUST include these files:

**Dockerfile** (multi-stage build):
```dockerfile
FROM node:20-alpine AS builder
WORKDIR /app
COPY package.json pnpm-lock.yaml ./
RUN corepack enable && pnpm install --frozen-lockfile
COPY . .
RUN pnpm build

FROM node:20-alpine
WORKDIR /app
COPY --from=builder /app/dist ./dist
COPY --from=builder /app/node_modules ./node_modules
COPY --from=builder /app/package.json ./
EXPOSE 3000
CMD ["node", "dist/index.js"]
```

**docker-compose.yml** (API + PostgreSQL + Redis):
```yaml
services:
  api:
    build: .
    ports:
      - "3000:3000"
    environment:
      - DATABASE_URL=postgres://postgres:postgres@db:5432/appdb?sslmode=disable
      - REDIS_URL=redis://redis:6379
      - JWT_SECRET=change-me
      - PORT=3000
    depends_on:
      db:
        condition: service_healthy
  db:
    image: postgres:16-alpine
    ports:
      - "5432:5432"
    environment:
      POSTGRES_DB: appdb
      POSTGRES_USER: postgres
      POSTGRES_PASSWORD: postgres
    volumes:
      - pgdata:/var/lib/postgresql/data
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U postgres"]
      interval: 5s
      timeout: 5s
      retries: 5
  redis:
    image: redis:7-alpine
    ports:
      - "6379:6379"
volumes:
  pgdata:
```

## Pre-Delivery Checklist

Before delivering Node.js/TypeScript code, verify:

- [ ] TypeScript strict mode enabled (tsconfig.json)
- [ ] No `any` types — use `unknown` + type guards
- [ ] All async errors handled (try/catch or .catch())
- [ ] Environment variables validated (Zod)
- [ ] Input validation on all endpoints (Zod)
- [ ] Graceful shutdown (SIGTERM/SIGINT)
- [ ] Structured logging (pino)
- [ ] Tests passing (vitest)
- [ ] ESLint + Prettier configured
- [ ] Dockerfile with multi-stage build
- [ ] No secrets in source code
- [ ] Database connections properly pooled
- [ ] Proper HTTP status codes used

## Quick Domain Reference

| Domain | When to Search |
|--------|---------------|
| `pattern` | Need a design pattern (factory, strategy, middleware, DI) |
| `package` | Choosing a library/package for a use case |
| `error` | Error handling strategy, custom errors, Result pattern |
| `perf` | Performance optimization, memory leaks, profiling |
| `test` | Testing approach, mocking, E2E, integration tests |
| `async` | Promises, streams, workers, event loop patterns |
| `typescript` | Type system patterns, generics, utility types |
| `arch` | Project structure, architecture decisions |
| `idiom` | Node.js/TypeScript conventions, naming, style |
| `anti` | Avoiding common mistakes and anti-patterns |
| `tooling` | Node.js tools, linting, profiling, CLI utilities |
| `module` | npm, packages, monorepo, ESM/CJS |

## Resources

- **scripts/**: BM25 search engine for intelligent knowledge base queries
- **data/**: 24 CSV databases with 590+ curated Node.js/TypeScript entries
- **Coverage**: 12 domains and 12 stacks for comprehensive Node.js development
