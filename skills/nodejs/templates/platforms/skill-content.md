# NodeSkill — Node.js/TypeScript Intelligence Toolkit

## Default Stack

| Component | Technology | Package |
|-----------|-----------|---------|
| Language | TypeScript 5.x | typescript |
| Runtime | Node.js 20+ | node |
| Web Framework | Express / Fastify | express, fastify |
| ORM | Prisma | @prisma/client |
| Validation | Zod | zod |
| Authentication | JWT | jsonwebtoken |
| Logging | Pino | pino, pino-http |
| Testing | Vitest | vitest |
| Deployment | Docker | Dockerfile |

## Search Commands

### Domain Search
```bash
# Design patterns (factory, strategy, middleware, etc.)
python3 {SKILL_DIR}/scripts/search.py "middleware pattern" --domain pattern

# Package recommendations
python3 {SKILL_DIR}/scripts/search.py "authentication library" --domain package

# Error handling strategies
python3 {SKILL_DIR}/scripts/search.py "custom error class" --domain error

# Performance optimization
python3 {SKILL_DIR}/scripts/search.py "event loop blocking" --domain perf

# Testing patterns
python3 {SKILL_DIR}/scripts/search.py "integration testing" --domain test

# Async/concurrency patterns
python3 {SKILL_DIR}/scripts/search.py "promise parallel" --domain async

# TypeScript type patterns
python3 {SKILL_DIR}/scripts/search.py "discriminated union" --domain typescript

# Architecture patterns
python3 {SKILL_DIR}/scripts/search.py "modular monolith" --domain arch

# Idioms and best practices
python3 {SKILL_DIR}/scripts/search.py "barrel exports" --domain idiom

# Anti-patterns to avoid
python3 {SKILL_DIR}/scripts/search.py "memory leak" --domain anti

# Tooling and CLI tools
python3 {SKILL_DIR}/scripts/search.py "eslint prettier" --domain tooling

# Module/package management
python3 {SKILL_DIR}/scripts/search.py "monorepo turborepo" --domain module
```

### Stack Search
```bash
# Express.js best practices
python3 {SKILL_DIR}/scripts/search.py "error middleware" --stack web-express

# Fastify patterns
python3 {SKILL_DIR}/scripts/search.py "plugin schema" --stack web-fastify

# NestJS patterns
python3 {SKILL_DIR}/scripts/search.py "guards interceptors" --stack web-nestjs

# Hono patterns
python3 {SKILL_DIR}/scripts/search.py "zod validator" --stack web-hono

# Next.js App Router
python3 {SKILL_DIR}/scripts/search.py "server components" --stack web-nextjs

# Remix patterns
python3 {SKILL_DIR}/scripts/search.py "loader action" --stack web-remix

# Prisma ORM
python3 {SKILL_DIR}/scripts/search.py "relations transactions" --stack orm-prisma

# Drizzle ORM
python3 {SKILL_DIR}/scripts/search.py "schema queries" --stack orm-drizzle

# GraphQL
python3 {SKILL_DIR}/scripts/search.py "dataloader resolver" --stack graphql

# Database patterns
python3 {SKILL_DIR}/scripts/search.py "connection pooling" --stack database

# Authentication
python3 {SKILL_DIR}/scripts/search.py "jwt refresh token" --stack auth

# Real-time communication
python3 {SKILL_DIR}/scripts/search.py "socket.io rooms" --stack realtime
```

### Auto-detect Domain
```bash
# No domain specified — auto-detects from query keywords
python3 {SKILL_DIR}/scripts/search.py "how to handle async errors in express"
```

## Pre-Delivery Checklist

Before delivering any Node.js/TypeScript code:
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
